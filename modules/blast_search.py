"""
BLAST Search Module
==================

Handles sequence similarity searches using BLAST and BLAT.
Includes local and remote BLAST searches, and result parsing.
"""

import sys
import os
import time
import logging
import requests
import tempfile
import subprocess
import re
from typing import Dict, List, Optional, Any
from urllib.parse import urlencode
from xml.etree import ElementTree as ET
from pathlib import Path

try:
    from Bio.Blast import NCBIWWW, NCBIXML
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    logging.warning("BioPython not available. Using basic sequence parsing.")
    # Minimal Seq and SeqRecord classes for fallback
    class Seq(str):
        def __new__(cls, sequence):
            return str.__new__(cls, sequence)
    class SeqRecord:
        def __init__(self, seq, id="", description="", letter_annotations=None):
            self.seq = seq
            self.id = id
            self.description = description
            self.letter_annotations = letter_annotations if letter_annotations is not None else {}
try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

class BlastSearcher:
    """
    BLAST and BLAT sequence similarity searcher.
    """
    
    def __init__(self, max_hits: int = 5, evalue_threshold: float = 1e-5, output_dir: str = None, force: bool = False):
        """
        Initialize BLAST searcher.
        
        Args:
            max_hits: Maximum number of hits to return per query
            evalue_threshold: E-value threshold for significance
            output_dir: Directory to save raw BLAST responses
            force: Force new BLAST searches even if cached results exist
        """
        self.max_hits = max_hits
        self.evalue_threshold = evalue_threshold
        self.output_dir = output_dir
        self.force = force
        self.logger = logging.getLogger(__name__)
        
        # Create raw_blast subdirectory if output_dir is provided
        if self.output_dir:
            self.raw_blast_dir = Path(output_dir) / "raw_blast_responses"
            self.raw_blast_dir.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Raw BLAST responses will be saved to: {self.raw_blast_dir}")
        else:
            self.raw_blast_dir = None
        self.logger.debug(f"BLAST searcher initialized with max_hits={self.max_hits}, evalue_threshold={self.evalue_threshold}, force={self.force}")
        # Check for local BLAST installation
        self.local_blast_available = self._check_local_blast()
        
    def _make_safe_filename(self, seq_id: str) -> str:
        self.logger.debug(f"Making safe filename for sequence ID: {seq_id}")
        """Convert sequence ID to safe filename."""
        # Replace problematic characters with underscores
        safe_name = re.sub(r'[<>:"/\\|?*]', '_', seq_id)
        # Limit length and remove any remaining problematic characters
        safe_name = safe_name[:50]  # Limit to 50 characters
        self.logger.debug(f"Safe filename: {safe_name}")
        return safe_name
        
    def _check_local_blast(self) -> bool:
        self.logger.debug("Checking for local BLAST+ installation...")
        """Check if local BLAST+ is available."""
        try:
            result = subprocess.run(['blastn', '-version'], 
                                  capture_output=True, text=True, timeout=10)
            self.logger.debug(f"Local BLAST+ version check output: {result.stdout}")
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            self.logger.debug("Local BLAST+ not found.")
            return False
    
    def blast_search(self, sequence: SeqRecord, database: str = 'nr', 
                    program: str = 'auto') -> Dict[str, Any]:
        """
        Perform BLAST search for a sequence.
        
        Args:
            sequence: SeqRecord object to search
            database: Database to search ('nr', 'nt', 'swissprot', etc.)
            program: BLAST program ('auto', 'blastn', 'blastp', 'blastx', etc.)
            
        Returns:
            Dictionary containing search results
        """
        self.logger.info(f"BLAST searching sequence: {sequence.id}")
        
        # Auto-detect BLAST program
        if program == 'auto':
            program = self._detect_blast_program(str(sequence.seq))

        # Check for cached results first (unless force is True)
        if not self.force and self.raw_blast_dir:
            self.logger.debug("Checking for cached BLAST results...")
            cached_result = self._check_cached_results(sequence.id, program, database)
            if cached_result:
                self.logger.debug(f"Using cached BLAST results for {sequence.id}")

                return cached_result
        
        try:
            return self._web_blast_search(sequence, database, program)
            # if BIOPYTHON_AVAILABLE:
            #     return self._biopython_blast_search(sequence, database, program)
        except Exception as e:
            self.logger.error(f"BLAST search failed: {str(e)}")
            return {'error': str(e), 'hits': []}
    
    def batch_blast_search(self, sequences: Dict[str, SeqRecord], database: str = 'nr', 
                          program: str = 'auto') -> Dict[str, Dict[str, Any]]:
        self.logger.debug("Starting batch BLAST search...")
        """
        Perform batch BLAST search for multiple sequences in a single submission.
        
        Args:
            sequences: Dictionary of sequence_id -> SeqRecord objects
            database: Database to search ('nr', 'nt', 'swissprot', etc.)
            program: BLAST program ('auto', 'blastn', 'blastp', 'blastx', etc.)
            
        Returns:
            Dictionary mapping sequence IDs to their search results
        """
        self.logger.info(f"Batch BLAST searching {len(sequences)} sequences...")
        
        # Auto-detect BLAST program from first sequence
        if program == 'auto':
            first_seq = next(iter(sequences.values()))
            program = self._detect_blast_program(str(first_seq.seq))
        
        # Check for cached results first
        results = {}
        uncached_sequences = {}
                        
        self.logger.debug("Checking for cached BLAST results in batch search...")
        if not self.force and self.raw_blast_dir:
            for seq_id, seq_record in sequences.items():
                cached_result = self._check_cached_results(seq_id, program, database)
                if cached_result:
                    self.logger.debug(f"Using cached BLAST results for {seq_id}")
                    results[seq_id] = cached_result
                else:
                    uncached_sequences[seq_id] = seq_record
        else:
            uncached_sequences = sequences
                                
        # If all results are cached, return them
        if not uncached_sequences:
            self.logger.info("All sequences have cached results")
            return results
        
        # Batch search uncached sequences
        self.logger.info(f"Performing batch BLAST search for {len(uncached_sequences)} sequences")
        
        try:
            batch_results = self._batch_web_blast_search(uncached_sequences, database, program)
            results.update(batch_results)
        except Exception as e:
            self.logger.error(f"Batch BLAST search failed, falling back to individual searches: {str(e)}")
            # Fallback to individual searches
            for seq_id, seq_record in uncached_sequences.items():
                results[seq_id] = self.blast_search(seq_record, database, program)
        
        return results
    
    def _detect_blast_program(self, sequence: str) -> str:
        """Auto-detect appropriate BLAST program based on sequence."""
        sequence = sequence.upper()
        nucleotide_chars = set('ATCGUN')
        seq_chars = set(sequence)
        
        # Check if mostly nucleotides
        nucleotide_fraction = len(seq_chars & nucleotide_chars) / len(seq_chars) if seq_chars else 0
        
        if nucleotide_fraction > 0.8:
            return 'blastn'  # Nucleotide to nucleotide
        else:
            return 'blastp'  # Protein to protein
    
    def _biopython_blast_search(self, sequence: SeqRecord, database: str, program: str) -> Dict[str, Any]:
        """Perform BLAST search using BioPython."""
        try:
            # Perform online BLAST search
            result_handle = NCBIWWW.qblast(
                program=program,
                database=database,
                sequence=str(sequence.seq),
                hitlist_size=self.max_hits,
                expect=self.evalue_threshold
            )
            
            # Save raw XML response if output directory is specified
            raw_xml_content = result_handle.read()
            if self.raw_blast_dir:
                safe_seq_id = self._make_safe_filename(sequence.id)
                raw_file = self.raw_blast_dir / f"{safe_seq_id}_biopython.xml"
                with open(raw_file, 'w') as f:
                    f.write(raw_xml_content)
                self.logger.info(f"Raw BioPython BLAST response saved: {raw_file}")
            
            # Parse results from the saved content
            from io import StringIO
            blast_records = NCBIXML.parse(StringIO(raw_xml_content))
            blast_record = next(blast_records)
            
            # Extract hits
            hits = []
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect <= self.evalue_threshold:
                        hit = {
                            'title': alignment.title,
                            'accession': alignment.accession,
                            'length': alignment.length,
                            'evalue': hsp.expect,
                            'score': hsp.score,
                            'bits': hsp.bits,
                            'identity': hsp.identities,
                            'positives': hsp.positives,
                            'gaps': hsp.gaps,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'subject_start': hsp.sbjct_start,
                            'subject_end': hsp.sbjct_end,
                            'query_sequence': hsp.query,
                            'subject_sequence': hsp.sbjct,
                            'match_sequence': hsp.match
                        }
                        hits.append(hit)
                        
                        if len(hits) >= self.max_hits:
                            break
                
                if len(hits) >= self.max_hits:
                    break
            
            return {
                'query_id': sequence.id,
                'program': program,
                'database': database,
                'hits': hits
            }
            
        except Exception as e:
            self.logger.error(f"BioPython BLAST search failed: {str(e)}")
            return {'error': str(e), 'hits': []}
    
    def _web_blast_search(self, sequence: SeqRecord, database: str, program: str) -> Dict[str, Any]:
        """Perform BLAST search using direct web interface (fallback)."""
        try:
            # Use NCBI BLAST web interface
            base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            # Submit search
            params = {
                'CMD': 'Put',
                'PROGRAM': program,
                'DATABASE': database,
                'QUERY': str(sequence.seq),
                'HITLIST_SIZE': self.max_hits,
                'EXPECT': self.evalue_threshold
            }
            
            response = requests.post(base_url, data=params, timeout=30)
            
            if response.status_code != 200:
                raise Exception(f"BLAST submission failed: {response.status_code}")
            
            # Extract RID (Request ID)
            rid = None
            for line in response.text.split('\n'):
                if 'RID = ' in line:
                    rid = line.split('RID = ')[1].strip()
                    break
            
            if not rid:
                raise Exception("Failed to get BLAST RID")
            
            # Poll for results
            self.logger.info(f"BLAST search submitted (RID: {rid}), waiting for results...")
            
            max_wait_time = 300  # 5 minutes
            start_time = time.time()
            
            while time.time() - start_time < max_wait_time:
                # Check status
                status_params = {'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid}
                status_response = requests.get(base_url, params=status_params, timeout=30)
                
                if 'Status=READY' in status_response.text:
                    break
                elif 'Status=WAITING' in status_response.text:
                    time.sleep(10)  # Wait 10 seconds before checking again
                else:
                    raise Exception("BLAST search failed or expired")
            
            else:
                raise Exception("BLAST search timed out")
            
            # Get results
            result_params = {
                'CMD': 'Get',
                'FORMAT_TYPE': 'XML',
                'RID': rid
            }
            
            result_response = requests.get(base_url, params=result_params, timeout=60)
            
            if result_response.status_code != 200:
                raise Exception(f"Failed to retrieve BLAST results: {result_response.status_code}")
            
            # Save raw XML response if output directory is specified
            raw_xml_content = result_response.text
            if self.raw_blast_dir:
                safe_seq_id = self._make_safe_filename(sequence.id)
                raw_file = self.raw_blast_dir / f"{safe_seq_id}_web.xml"
                with open(raw_file, 'w') as f:
                    f.write(raw_xml_content)
                self.logger.info(f"Raw web BLAST response saved: {raw_file}")
            
            # Parse XML results
            hits = self._parse_blast_xml(raw_xml_content)
            
            return {
                'query_id': sequence.id,
                'program': program,
                'database': database,
                'hits': hits,
                'rid': rid
            }
            
        except Exception as e:
            self.logger.error(f"Web BLAST search failed: {str(e)}")
            return {'error': str(e), 'hits': []}
    
    def _batch_web_blast_search(self, sequences: Dict[str, SeqRecord], database: str, program: str) -> Dict[str, Dict[str, Any]]:
        """Perform batch BLAST search using direct web interface."""
        try:
            # Use NCBI BLAST web interface
            base_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
            
            # Create multi-FASTA query string
            query_fasta = ""
            seq_order = []  # Keep track of sequence order
            
            for seq_id, seq_record in sequences.items():
                query_fasta += f">{seq_id}\n{str(seq_record.seq)}\n"
                seq_order.append(seq_id)
            
            self.logger.debug(f"Created multi-FASTA query with {len(sequences)} sequences")
            
            # Submit batch search
            params = {
                'CMD': 'Put',
                'PROGRAM': program,
                'DATABASE': database,
                'QUERY': query_fasta,
                'HITLIST_SIZE': self.max_hits,
                'EXPECT': self.evalue_threshold
            }
            
            response = requests.post(base_url, data=params, timeout=60)
            
            if response.status_code != 200:
                raise Exception(f"Batch BLAST submission failed: {response.status_code}")
            
            # Extract RID (Request ID)
            rid = None
            for line in response.text.split('\n'):
                if 'RID = ' in line:
                    rid = line.split('RID = ')[1].strip()
                    break
            
            if not rid:
                raise Exception("Failed to get batch BLAST RID")
            
            # Poll for results
            self.logger.info(f"Batch BLAST search submitted (RID: {rid}), waiting for results...")
            
            max_wait_time = 600  # 10 minutes for batch jobs
            start_time = time.time()
            
            while time.time() - start_time < max_wait_time:
                # Check status
                status_params = {'CMD': 'Get', 'FORMAT_OBJECT': 'SearchInfo', 'RID': rid}
                status_response = requests.get(base_url, params=status_params, timeout=30)
                
                if 'Status=READY' in status_response.text:
                    break
                elif 'Status=WAITING' in status_response.text:
                    time.sleep(15)  # Wait longer for batch jobs
                else:
                    raise Exception("Batch BLAST search failed or expired")
            
            else:
                raise Exception("Batch BLAST search timed out")
            
            # Get results
            result_params = {
                'CMD': 'Get',
                'FORMAT_TYPE': 'XML',
                'RID': rid
            }
            
            result_response = requests.get(base_url, params=result_params, timeout=120)
            
            if result_response.status_code != 200:
                raise Exception(f"Failed to retrieve batch BLAST results: {result_response.status_code}")
            
            # Save raw XML response if output directory is specified
            raw_xml_content = result_response.text
            if self.raw_blast_dir:
                raw_file = self.raw_blast_dir / f"batch_{rid}_web.xml"
                with open(raw_file, 'w') as f:
                    f.write(raw_xml_content)
                self.logger.info(f"Raw batch BLAST response saved: {raw_file}")
            
            # Parse batch XML results
            batch_results = self._parse_batch_blast_xml(raw_xml_content, seq_order, program, database, rid)
            
            return batch_results
            
        except Exception as e:
            self.logger.error(f"Batch web BLAST search failed: {str(e)}")
            raise e
    
    def _parse_blast_xml(self, xml_content: str) -> List[Dict[str, Any]]:
        """Parse BLAST XML results."""
        hits = []
        
        try:
            root = ET.fromstring(xml_content)
            
            # Find iterations - use direct path for better reliability
            iterations = root.findall('.//Iteration')
            self.logger.debug(f"Found {len(iterations)} iterations in XML")
            
            for iteration in iterations:
                hit_elements = iteration.findall('.//Hit')
                self.logger.debug(f"Found {len(hit_elements)} hits in iteration")
                
                for hit in hit_elements:
                    hit_def = hit.find('Hit_def')
                    hit_accession = hit.find('Hit_accession')
                    hit_len = hit.find('Hit_len')
                    
                    # Process HSPs (High-scoring Segment Pairs)
                    hsps = hit.findall('.//Hsp')
                    for hsp in hsps:
                        try:
                            evalue_elem = hsp.find('Hsp_evalue')
                            if evalue_elem is None:
                                continue
                                
                            evalue = float(evalue_elem.text)
                            self.logger.debug(f"Processing hit with e-value: {evalue}, threshold: {self.evalue_threshold}")
                            
                            # Apply e-value filter but also log what we're finding
                            if evalue <= self.evalue_threshold:
                                hit_data = {
                                    'title': hit_def.text if hit_def is not None else 'Unknown',
                                    'accession': hit_accession.text if hit_accession is not None else 'Unknown',
                                    'length': int(hit_len.text) if hit_len is not None else 0,
                                    'evalue': evalue,
                                    'score': float(hsp.find('Hsp_score').text),
                                    'bits': float(hsp.find('Hsp_bit-score').text),
                                    'identity': int(hsp.find('Hsp_identity').text),
                                    'positives': int(hsp.find('Hsp_positive').text),
                                    'gaps': int(hsp.find('Hsp_gaps').text),
                                    'query_start': int(hsp.find('Hsp_query-from').text),
                                    'query_end': int(hsp.find('Hsp_query-to').text),
                                    'subject_start': int(hsp.find('Hsp_hit-from').text),
                                    'subject_end': int(hsp.find('Hsp_hit-to').text),
                                    'query_sequence': hsp.find('Hsp_qseq').text,
                                    'subject_sequence': hsp.find('Hsp_hseq').text,
                                    'match_sequence': hsp.find('Hsp_midline').text
                                }
                                
                                hits.append(hit_data)
                                self.logger.debug(f"Added hit: {hit_data['title'][:50]}... (e-value: {evalue})")
                                
                                if len(hits) >= self.max_hits:
                                    break
                            else:
                                self.logger.debug(f"Skipping hit due to e-value {evalue} > {self.evalue_threshold}")
                        except (ValueError, AttributeError) as e:
                            self.logger.warning(f"Error parsing HSP: {e}")
                            continue
                    
                    if len(hits) >= self.max_hits:
                        break
        
        except ET.ParseError as e:
            self.logger.error(f"Failed to parse BLAST XML: {str(e)}")
        except Exception as e:
            self.logger.error(f"Unexpected error parsing BLAST XML: {str(e)}")
        
        self.logger.info(f"Parsed {len(hits)} valid hits from XML")
        return hits
    
    def _parse_batch_blast_xml(self, xml_content: str, seq_order: List[str], 
                              program: str, database: str, rid: str) -> Dict[str, Dict[str, Any]]:
        """Parse batch BLAST XML results and separate by query sequence."""
        results = {}
        
        try:
            root = ET.fromstring(xml_content)
            
            # Find all iterations - each iteration corresponds to one query sequence
            iterations = root.findall('.//Iteration')
            self.logger.debug(f"Found {len(iterations)} iterations in batch XML")
            
            for i, iteration in enumerate(iterations):
                # Get query information
                iter_query = iteration.find('Iteration_query-def')
                query_id = iter_query.text if iter_query is not None else f"query_{i}"
                
                # Sometimes BLAST modifies the query ID, so use order if possible
                if i < len(seq_order):
                    query_id = seq_order[i]
                
                self.logger.debug(f"Processing iteration {i}: {query_id}")
                
                # Parse hits for this iteration
                hits = []
                hit_elements = iteration.findall('.//Hit')
                self.logger.debug(f"Found {len(hit_elements)} hits for {query_id}")
                
                for hit in hit_elements:
                    hit_def = hit.find('Hit_def')
                    hit_accession = hit.find('Hit_accession')
                    hit_len = hit.find('Hit_len')
                    
                    # Process HSPs (High-scoring Segment Pairs)
                    hsps = hit.findall('.//Hsp')
                    for hsp in hsps:
                        try:
                            evalue_elem = hsp.find('Hsp_evalue')
                            if evalue_elem is None:
                                continue
                                
                            evalue = float(evalue_elem.text)
                            
                            # Apply e-value filter
                            if evalue <= self.evalue_threshold:
                                hit_data = {
                                    'title': hit_def.text if hit_def is not None else 'Unknown',
                                    'accession': hit_accession.text if hit_accession is not None else 'Unknown',
                                    'length': int(hit_len.text) if hit_len is not None else 0,
                                    'evalue': evalue,
                                    'score': float(hsp.find('Hsp_score').text),
                                    'bits': float(hsp.find('Hsp_bit-score').text),
                                    'identity': int(hsp.find('Hsp_identity').text),
                                    'positives': int(hsp.find('Hsp_positive').text),
                                    'gaps': int(hsp.find('Hsp_gaps').text),
                                    'query_start': int(hsp.find('Hsp_query-from').text),
                                    'query_end': int(hsp.find('Hsp_query-to').text),
                                    'subject_start': int(hsp.find('Hsp_hit-from').text),
                                    'subject_end': int(hsp.find('Hsp_hit-to').text),
                                    'query_sequence': hsp.find('Hsp_qseq').text,
                                    'subject_sequence': hsp.find('Hsp_hseq').text,
                                    'match_sequence': hsp.find('Hsp_midline').text
                                }
                                
                                hits.append(hit_data)
                                
                                if len(hits) >= self.max_hits:
                                    break
                        except (ValueError, AttributeError) as e:
                            self.logger.warning(f"Error parsing HSP for {query_id}: {e}")
                            continue
                    
                    if len(hits) >= self.max_hits:
                        break
                
                # Store results for this query
                results[query_id] = {
                    'query_id': query_id,
                    'program': program,
                    'database': database,
                    'hits': hits,
                    'rid': rid,
                    'batch': True
                }
                
                self.logger.debug(f"Parsed {len(hits)} hits for {query_id}")
                
                # Save individual results to cache if requested
                if self.raw_blast_dir and hits:
                    safe_seq_id = self._make_safe_filename(query_id)
                    cache_file = self.raw_blast_dir / f"{safe_seq_id}_batch.xml"
                    
                    # Create individual XML for caching
                    individual_xml = self._create_individual_xml_from_iteration(iteration, query_id, rid)
                    with open(cache_file, 'w') as f:
                        f.write(individual_xml)
                    self.logger.debug(f"Cached individual result: {cache_file}")
        
        except ET.ParseError as e:
            self.logger.error(f"Failed to parse batch BLAST XML: {str(e)}")
        except Exception as e:
            self.logger.error(f"Unexpected error parsing batch BLAST XML: {str(e)}")
        
        self.logger.info(f"Parsed batch results for {len(results)} sequences")
        return results
    
    def _create_individual_xml_from_iteration(self, iteration, query_id: str, rid: str) -> str:
        """Create individual XML file content from batch iteration for caching."""
        # Create a minimal XML structure that matches single-query format
        xml_template = f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLAST 2.0+</BlastOutput_version>
  <BlastOutput_reference>Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.</BlastOutput_reference>
  <BlastOutput_db>nr</BlastOutput_db>
  <BlastOutput_query-ID>{query_id}</BlastOutput_query-ID>
  <BlastOutput_query-def>{query_id}</BlastOutput_query-def>
  <BlastOutput_query-len>0</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_expect>{self.evalue_threshold}</Parameters_expect>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    {ET.tostring(iteration, encoding='unicode')}
  </BlastOutput_iterations>
</BlastOutput>"""
        return xml_template
    
    def _check_cached_results(self, sequence_id: str, program: str, database: str) -> Optional[Dict[str, Any]]:
        """
        Check for cached BLAST results and parse them if found.
        
        Args:
            sequence_id: ID of the sequence
            program: BLAST program used
            database: Database searched
            
        Returns:
            Parsed results dictionary if cached results exist, None otherwise
        """
        if not self.raw_blast_dir:
            return None
            
        safe_seq_id = self._make_safe_filename(sequence_id)
        
        # Check for both BioPython and web BLAST cached files
        cached_files = [
            self.raw_blast_dir / f"{safe_seq_id}_biopython.xml",
            self.raw_blast_dir / f"{safe_seq_id}_web.xml",
            self.raw_blast_dir / f"{safe_seq_id}_batch.xml"
        ]
        self.logger.debug(f"Checking for cached files: {cached_files}")
        for cached_file in cached_files:
            self.logger.debug(f"Looking for cached file: {cached_file}")
            if cached_file.exists():
                try:
                    with open(cached_file, 'r') as f:
                        xml_content = f.read()
                    self.logger.debug(f"Found cached file: {cached_file}")
                    # Parse the cached XML
                    hits = self._parse_blast_xml(xml_content)
                    
                    # Determine which method was used based on filename
                    method = "biopython" if "_biopython.xml" in str(cached_file) else "web" if "_web.xml" in str(cached_file) else "batch" if "_batch.xml" in str(cached_file) else "unknown"
                    self.logger.debug(f"Parsed {len(hits)} hits from cached file: {cached_file}")
                    return {
                        'query_id': sequence_id,
                        'program': program,
                        'database': database,
                        'hits': hits,
                        'cached': True,
                        'cache_file': str(cached_file),
                        'method': method
                    }
                    
                except Exception as e:
                    self.logger.warning(f"Failed to parse cached results from {cached_file}: {e}")
                    continue
        self.logger.debug(f"No cached results found for {sequence_id}")
        return None
    
    def extract_organism_info(self, hit: Dict[str, Any]) -> Dict[str, str]:
        """
        Extract organism information from BLAST hit.
        
        Args:
            hit: BLAST hit dictionary
            
        Returns:
            Dictionary with organism information
        """
        title = hit.get('title', '')
        
        # Common patterns for organism extraction
        organism_info = {
            'organism': 'Unknown',
            'description': title,
            'gene': 'Unknown'
        }
        
        # Try to extract organism name
        import re
        
        # Pattern 1: [Organism name]
        bracket_match = re.search(r'\[([^\]]+)\]', title)
        if bracket_match:
            organism_info['organism'] = bracket_match.group(1)
        
        # Pattern 2: Common prefixes
        for pattern in ['OS=', 'organism=']:
            if pattern in title.lower():
                parts = title.lower().split(pattern)
                if len(parts) > 1:
                    organism_part = parts[1].split()[0:2]  # Take first two words
                    organism_info['organism'] = ' '.join(organism_part)
                    break
        
        # Try to extract gene information
        gene_patterns = ['gene=', 'GN=']
        for pattern in gene_patterns:
            if pattern in title:
                parts = title.split(pattern)
                if len(parts) > 1:
                    gene_part = parts[1].split()[0]
                    organism_info['gene'] = gene_part
                    break
        
        return organism_info
    
    def summarize_search_results(self, search_results: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Summarize multiple search results.
        
        Args:
            search_results: List of search result dictionaries
            
        Returns:
            Summary dictionary
        """
        if not search_results:
            return {'summary': 'No search results to summarize'}
        
        # Collect all hits
        all_hits = []
        successful_searches = 0
        
        for result in search_results:
            if 'hits' in result and result['hits']:
                all_hits.extend(result['hits'])
                successful_searches += 1
        
        if not all_hits:
            return {'summary': 'No significant hits found in any search'}
        
        # Extract organism information
        organisms = []
        for hit in all_hits:
            org_info = self.extract_organism_info(hit)
            organisms.append(org_info['organism'])
        
        # Count organism occurrences
        organism_counts = {}
        for org in organisms:
            organism_counts[org] = organism_counts.get(org, 0) + 1
        
        # Sort by frequency
        top_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        # Calculate statistics
        evalues = [hit['evalue'] for hit in all_hits if 'evalue' in hit]
        identities = [hit.get('identity', 0) for hit in all_hits]
        
        summary = {
            'total_searches': len(search_results),
            'successful_searches': successful_searches,
            'total_hits': len(all_hits),
            'top_organisms': top_organisms,
            'min_evalue': min(evalues) if evalues else None,
            'max_evalue': max(evalues) if evalues else None,
            'mean_identity': sum(identities) / len(identities) if identities else 0,
            'max_identity': max(identities) if identities else 0
        }
        
        return summary
    
    def create_results_dataframe(self, search_results: List[Dict[str, Any]]):
        """
        Create pandas DataFrame from search results.
        
        Args:
            search_results: List of search result dictionaries
            
        Returns:
            pandas DataFrame or None if pandas not available
        """
        if not PANDAS_AVAILABLE:
            self.logger.warning("Pandas not available. Cannot create DataFrame.")
            return None
        
        # Flatten results into rows
        rows = []
        
        for result in search_results:
            query_id = result.get('query_id', 'Unknown')
            
            if 'hits' in result and result['hits']:
                for i, hit in enumerate(result['hits'][:self.max_hits]):
                    org_info = self.extract_organism_info(hit)
                    
                    row = {
                        'query_id': query_id,
                        'hit_rank': i + 1,
                        'title': hit.get('title', ''),
                        'accession': hit.get('accession', ''),
                        'organism': org_info['organism'],
                        'gene': org_info['gene'],
                        'evalue': hit.get('evalue', None),
                        'identity': hit.get('identity', 0),
                        'score': hit.get('score', 0),
                        'length': hit.get('length', 0)
                    }
                    rows.append(row)
            else:
                # No hits found
                row = {
                    'query_id': query_id,
                    'hit_rank': 0,
                    'title': 'No significant hits',
                    'accession': '',
                    'organism': 'Unknown',
                    'gene': 'Unknown',
                    'evalue': None,
                    'identity': 0,
                    'score': 0,
                    'length': 0
                }
                rows.append(row)
        
        df = pd.DataFrame(rows)
        return df