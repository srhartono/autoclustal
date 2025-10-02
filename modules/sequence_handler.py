"""
Sequence Handler Module
======================

Handles loading, parsing, and preprocessing of FASTA/FASTQ sequence files.
Includes automatic sequence type detection and validation.
"""

import os
import re
import glob
import logging
import random
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from collections import defaultdict

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
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

class SequenceHandler:
    """
    Handler for sequence file I/O and preprocessing.
    """
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.supported_formats = ['fasta', 'fastq', 'fa', 'fq', 'fas']
        
    def load_sequences(self, input_files: List[str], seq_type: str = 'auto', 
                      max_sequences: Optional[int] = None, 
                      random_sample_size: Optional[int] = None,
                      random_seed: Optional[int] = None) -> Dict[str, SeqRecord]:
        """
        Load sequences from input files with optional sampling.
        
        Args:
            input_files: List of file paths (supports wildcards)
            seq_type: Sequence type ('auto', 'dna', 'rna', 'protein')
            max_sequences: Maximum number of sequences to take from the beginning (M parameter)
            random_sample_size: Number of random sequences to sample from the first max_sequences (N parameter)
            random_seed: Random seed for reproducible sampling
            
        Returns:
            Dictionary of sequence ID to SeqRecord objects
        """
        sequences = {}
        file_paths = self._expand_file_paths(input_files)
        
        # Set random seed for reproducibility
        if random_seed is not None:
            random.seed(random_seed)
        
        all_sequences = []
        
        for file_path in file_paths:
            self.logger.info(f"Loading sequences from: {file_path}")
            file_sequences = self._load_single_file(file_path)
            all_sequences.extend(file_sequences)
        
        # Apply sampling if specified
        if max_sequences is not None:
            if len(all_sequences) > max_sequences:
                self.logger.info(f"Taking first {max_sequences} sequences from {len(all_sequences)} total")
                all_sequences = all_sequences[:max_sequences]
        
        if random_sample_size is not None:
            if len(all_sequences) > random_sample_size:
                self.logger.info(f"Randomly sampling {random_sample_size} sequences from {len(all_sequences)} available")
                all_sequences = random.sample(all_sequences, random_sample_size)
            else:
                self.logger.info(f"Requested {random_sample_size} random samples, but only {len(all_sequences)} sequences available. Using all sequences.")
        
        # Detect sequence type if auto
        if seq_type == 'auto':
            detected_type = self._detect_sequence_type(all_sequences)
            self.logger.info(f"Detected sequence type: {detected_type}")
        else:
            detected_type = seq_type
        
        # Validate and process sequences
        processed_sequences = self._process_sequences(all_sequences, detected_type)
        sequences.update(processed_sequences)
        
        self.logger.info(f"Loaded {len(sequences)} total sequences after sampling")
        return sequences
    
    def _expand_file_paths(self, input_files: List[str]) -> List[str]:
        """Expand wildcards in file paths."""
        file_paths = []
        for pattern in input_files:
            if '*' in pattern or '?' in pattern:
                expanded = glob.glob(pattern)
                file_paths.extend(expanded)
            else:
                file_paths.append(pattern)
        
        # Filter existing files
        existing_files = [f for f in file_paths if os.path.exists(f)]
        if not existing_files:
            raise FileNotFoundError(f"No valid files found: {input_files}")
        
        return existing_files
    
    def _load_single_file(self, file_path: str) -> List[SeqRecord]:
        """Load sequences from a single file."""
        file_format = self._detect_file_format(file_path)
        sequences = []
        
        if BIOPYTHON_AVAILABLE:
            try:
                sequences = list(SeqIO.parse(file_path, file_format))
            except Exception as e:
                self.logger.warning(f"BioPython parsing failed: {e}. Using basic parser.")
                sequences = self._basic_parse(file_path, file_format)
        else:
            sequences = self._basic_parse(file_path, file_format)
        
        return sequences
    
    def _detect_file_format(self, file_path: str) -> str:
        """Detect file format from extension and content."""
        extension = Path(file_path).suffix.lower().lstrip('.')
        
        if extension in ['fastq', 'fq']:
            return 'fastq'
        elif extension in ['fasta', 'fa', 'fas', 'fna', 'faa']:
            return 'fasta'
        else:
            # Try to detect from content
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                if first_line.startswith('@'):
                    return 'fastq'
                elif first_line.startswith('>'):
                    return 'fasta'
        
        # Default to fasta
        return 'fasta'
    
    def _basic_parse(self, file_path: str, file_format: str) -> List[SeqRecord]:
        """Basic sequence parser when BioPython is not available."""
        sequences = []
        
        with open(file_path, 'r') as f:
            if file_format == 'fasta':
                sequences = self._parse_fasta(f)
            elif file_format == 'fastq':
                sequences = self._parse_fastq(f)
        
        return sequences
    
    def _parse_fasta(self, file_handle) -> List[SeqRecord]:
        """Parse FASTA format."""
        sequences = []
        current_id = None
        current_seq = []
        
        for line in file_handle:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    seq_record = SeqRecord(
                        Seq(''.join(current_seq)),
                        id=current_id,
                        description=''
                    )
                    sequences.append(seq_record)
                
                current_id = line[1:].split()[0]  # Take first word after >
                current_seq = []
            else:
                current_seq.append(line)
        
        # Add last sequence
        if current_id:
            seq_record = SeqRecord(
                Seq(''.join(current_seq)),
                id=current_id,
                description=''
            )
            sequences.append(seq_record)
        
        return sequences
    
    def _parse_fastq(self, file_handle) -> List[SeqRecord]:
        """Parse FASTQ format."""
        sequences = []
        lines = file_handle.readlines()
        
        for i in range(0, len(lines), 4):
            if i + 3 < len(lines):
                header = lines[i].strip()
                sequence = lines[i + 1].strip()
                quality = lines[i + 3].strip()
                
                seq_id = header[1:].split()[0]  # Remove @ and take first word
                seq_record = SeqRecord(
                    Seq(sequence),
                    id=seq_id,
                    description='',
                    letter_annotations={'phred_quality': [ord(c) - 33 for c in quality]}
                )
                sequences.append(seq_record)
        
        return sequences
    
    def _detect_sequence_type(self, sequences: List[SeqRecord]) -> str:
        """
        Detect sequence type based on nucleotide composition.
        
        Args:
            sequences: List of SeqRecord objects
            
        Returns:
            Detected sequence type: 'dna', 'rna', or 'protein'
        """
        if not sequences:
            return 'dna'  # Default
        
        # Sample up to 100 sequences for type detection
        sample_sequences = sequences[:min(100, len(sequences))]
        total_length = 0
        nucleotide_counts = defaultdict(int)
        
        for seq_record in sample_sequences:
            sequence = str(seq_record.seq).upper()
            total_length += len(sequence)
            
            for nucleotide in sequence:
                nucleotide_counts[nucleotide] += 1
        
        # Calculate composition
        if total_length == 0:
            return 'dna'
        
        # Check for protein-specific amino acids
        protein_specific = set('EFILPQZ')
        dna_nucleotides = set('ATCG')
        rna_nucleotides = set('AUCG')
        
        # Check if we have protein-specific amino acids
        if any(aa in nucleotide_counts for aa in protein_specific):
            return 'protein'
        
        # Check nucleotide composition
        total_nucleotides = sum(nucleotide_counts[nt] for nt in 'ATUCGN')
        nucleotide_fraction = total_nucleotides / total_length
        
        if nucleotide_fraction > 0.9:  # Mostly nucleotides
            # Check for RNA (presence of U)
            if nucleotide_counts['U'] > nucleotide_counts['T']:
                return 'rna'
            else:
                return 'dna'
        else:
            return 'protein'
    
    def _process_sequences(self, sequences: List[SeqRecord], seq_type: str) -> Dict[str, SeqRecord]:
        """
        Process and validate sequences.
        
        Args:
            sequences: List of SeqRecord objects
            seq_type: Sequence type
            
        Returns:
            Dictionary of processed sequences
        """
        processed = {}
        
        for i, seq_record in enumerate(sequences):
            # Ensure unique IDs
            seq_id = seq_record.id if seq_record.id else f"seq_{i+1}"
            if seq_id in processed:
                seq_id = f"{seq_id}_{i+1}"
            
            # Clean sequence
            sequence = str(seq_record.seq).upper()
            sequence = re.sub(r'[^A-Z]', '', sequence)  # Remove non-alphabetic characters
            
            # Validate sequence
            if self._validate_sequence(sequence, seq_type):
                processed_record = SeqRecord(
                    Seq(sequence),
                    id=seq_id,
                    description=seq_record.description
                )
                processed[seq_id] = processed_record
            else:
                self.logger.warning(f"Skipping invalid sequence: {seq_id}")
        
        return processed
    
    def _validate_sequence(self, sequence: str, seq_type: str) -> bool:
        """
        Validate sequence based on type.
        
        Args:
            sequence: Sequence string
            seq_type: Expected sequence type
            
        Returns:
            True if valid, False otherwise
        """
        if len(sequence) < 10:  # Minimum length requirement
            return False
        
        valid_chars = {
            'dna': set('ATCGN'),
            'rna': set('AUCGN'),
            'protein': set('ACDEFGHIKLMNPQRSTVWY')
        }
        
        if seq_type in valid_chars:
            sequence_chars = set(sequence)
            invalid_chars = sequence_chars - valid_chars[seq_type]
            
            # Allow up to 5% invalid characters (for ambiguous bases/residues)
            invalid_fraction = len([c for c in sequence if c in invalid_chars]) / len(sequence)
            return invalid_fraction <= 0.05
        
        return True
    
    def save_sequences(self, sequences: Dict[str, SeqRecord], output_path: str, file_format: str = 'fasta'):
        """
        Save sequences to file.
        
        Args:
            sequences: Dictionary of sequences
            output_path: Output file path
            file_format: Output format ('fasta' or 'fastq')
        """
        sequence_list = list(sequences.values())
        
        if BIOPYTHON_AVAILABLE:
            SeqIO.write(sequence_list, output_path, file_format)
        else:
            self._basic_write(sequence_list, output_path, file_format)
        
        self.logger.info(f"Saved {len(sequences)} sequences to {output_path}")
    
    def _basic_write(self, sequences: List[SeqRecord], output_path: str, file_format: str):
        """Basic sequence writer when BioPython is not available."""
        with open(output_path, 'w') as f:
            for seq_record in sequences:
                if file_format == 'fasta':
                    f.write(f">{seq_record.id}\n")
                    # Write sequence in lines of 80 characters
                    sequence = str(seq_record.seq)
                    for i in range(0, len(sequence), 80):
                        f.write(f"{sequence[i:i+80]}\n")
                elif file_format == 'fastq':
                    f.write(f"@{seq_record.id}\n")
                    f.write(f"{seq_record.seq}\n")
                    f.write("+\n")
                    # Use default quality scores if not available
                    quality = getattr(seq_record, 'letter_annotations', {}).get('phred_quality', [40] * len(seq_record.seq))
                    quality_string = ''.join(chr(q + 33) for q in quality)
                    f.write(f"{quality_string}\n")
    
    def get_sequence_stats(self, sequences: Dict[str, SeqRecord]) -> Dict[str, float]:
        """
        Calculate basic sequence statistics.
        
        Args:
            sequences: Dictionary of sequences
            
        Returns:
            Dictionary of statistics
        """
        if not sequences:
            return {}
        
        lengths = [len(seq.seq) for seq in sequences.values()]
        
        stats = {
            'count': len(sequences),
            'min_length': min(lengths),
            'max_length': max(lengths),
            'mean_length': sum(lengths) / len(lengths),
            'total_length': sum(lengths)
        }
        
        return stats
    
    def sample_sequences_from_top(self, sequences: List[SeqRecord], 
                                 max_sequences: int, 
                                 random_sample_size: int,
                                 random_seed: Optional[int] = None) -> List[SeqRecord]:
        """
        Sample random N sequences from top M sequences.
        
        Args:
            sequences: List of SeqRecord objects
            max_sequences: Maximum number of sequences to consider from the top (M)
            random_sample_size: Number of sequences to randomly sample (N)
            random_seed: Random seed for reproducibility
            
        Returns:
            List of sampled SeqRecord objects
        """
        if random_seed is not None:
            random.seed(random_seed)
        
        # Take first M sequences
        top_sequences = sequences[:max_sequences] if len(sequences) > max_sequences else sequences
        
        # Sample N sequences randomly from the top M
        if len(top_sequences) <= random_sample_size:
            self.logger.warning(f"Requested {random_sample_size} samples from {len(top_sequences)} sequences. Returning all available sequences.")
            return top_sequences
        
        sampled = random.sample(top_sequences, random_sample_size)
        self.logger.info(f"Sampled {len(sampled)} sequences from top {len(top_sequences)} sequences")
        
        return sampled