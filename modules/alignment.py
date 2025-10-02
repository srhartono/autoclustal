"""
Alignment Module
===============

Handles multiple sequence alignment using various alignment tools.
Supports MUSCLE, ClustalW, and MAFFT aligners.
"""

import os
import subprocess
import tempfile
import logging
from pathlib import Path
from typing import Dict, List, Optional

try:
    from Bio import SeqIO, AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment
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
    class MultipleSeqAlignment:
        def __init__(self, records):
            self._records = records
        def __iter__(self):
            return iter(self._records)
        def __len__(self):
            return len(self._records)
        def __getitem__(self, index):
            return self._records[index]

class AlignmentEngine:
    """
    Multiple sequence alignment engine supporting various alignment tools.
    """
    
    def __init__(self, tool: str = 'muscle', threads: int = 4):
        """
        Initialize alignment engine.
        
        Args:
            tool: Alignment tool ('muscle', 'clustalw', 'mafft')
            threads: Number of threads to use
        """
        self.tool = tool.lower()
        self.threads = threads
        self.logger = logging.getLogger(__name__)
        
        # Check if alignment tool is available
        self.tool_available = self._check_tool_availability()
        
        if not self.tool_available:
            self.logger.warning(f"{self.tool} not found. Will use basic pairwise alignment.")
    
    def _check_tool_availability(self) -> bool:
        """Check if the selected alignment tool is available."""
        tool_commands = {
            'muscle': 'muscle',
            'clustalw': 'clustalw',
            'mafft': 'mafft'
        }
        
        command = tool_commands.get(self.tool)
        if not command:
            return False
        
        try:
            result = subprocess.run([command, '--help'], 
                                  capture_output=True, 
                                  text=True, 
                                  timeout=10)
            return result.returncode == 0 or 'usage' in result.stderr.lower()
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def align_sequences(self, sequences: Dict[str, SeqRecord], output_path: Optional[str] = None) -> MultipleSeqAlignment:
        """
        Perform multiple sequence alignment.
        
        Args:
            sequences: Dictionary of sequence ID to SeqRecord
            output_path: Output file path for alignment
            
        Returns:
            MultipleSeqAlignment object or list of SeqRecords
        """
        self.logger.info(f"Aligning {len(sequences)} sequences using {self.tool}")
        
        if len(sequences) < 2:
            raise ValueError("At least 2 sequences required for alignment")
        
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
            temp_input_path = temp_input.name
            self._write_sequences_to_file(sequences, temp_input_path)
        
        try:
            if self.tool_available:
                alignment = self._run_external_aligner(temp_input_path)
            else:
                # Fallback to basic alignment
                alignment = self._basic_alignment(sequences)
            
            # Save alignment if output path provided
            if output_path:
                self._save_alignment(alignment, output_path)
            
            return alignment
            
        finally:
            # Clean up temporary file
            if os.path.exists(temp_input_path):
                os.unlink(temp_input_path)
    
    def _write_sequences_to_file(self, sequences: Dict[str, SeqRecord], file_path: str):
        """Write sequences to FASTA file."""
        with open(file_path, 'w') as f:
            for seq_record in sequences.values():
                f.write(f">{seq_record.id}\n")
                # Write sequence in lines of 80 characters
                sequence = str(seq_record.seq)
                for i in range(0, len(sequence), 80):
                    f.write(f"{sequence[i:i+80]}\n")
    
    def _run_external_aligner(self, input_path: str) -> MultipleSeqAlignment:
        """Run external alignment tool."""
        with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        try:
            if self.tool == 'muscle':
                self._run_muscle(input_path, temp_output_path)
            elif self.tool == 'clustalw':
                self._run_clustalw(input_path, temp_output_path)
            elif self.tool == 'mafft':
                self._run_mafft(input_path, temp_output_path)
            
            # Read alignment result
            if BIOPYTHON_AVAILABLE:
                alignment = AlignIO.read(temp_output_path, 'fasta')
            else:
                alignment = self._read_alignment_basic(temp_output_path)
            
            return alignment
            
        finally:
            if os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
    
    def _run_muscle(self, input_path: str, output_path: str):
        """Run MUSCLE alignment."""
        cmd = ['muscle', '-in', input_path, '-out', output_path]
        if self.threads > 1:
            cmd.extend(['-threads', str(self.threads)])
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MUSCLE failed: {result.stderr}")
    
    def _run_clustalw(self, input_path: str, output_path: str):
        """Run ClustalW alignment."""
        # ClustalW creates output with specific naming
        base_name = os.path.splitext(input_path)[0]
        clustal_output = f"{base_name}.aln"
        
        cmd = ['clustalw', '-infile=' + input_path, '-outfile=' + clustal_output]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"ClustalW failed: {result.stderr}")
        
        # Convert to FASTA format
        if os.path.exists(clustal_output):
            if BIOPYTHON_AVAILABLE:
                alignment = AlignIO.read(clustal_output, 'clustal')
                AlignIO.write(alignment, output_path, 'fasta')
            else:
                # Basic conversion from clustal to fasta
                self._convert_clustal_to_fasta(clustal_output, output_path)
            os.unlink(clustal_output)
    
    def _run_mafft(self, input_path: str, output_path: str):
        """Run MAFFT alignment."""
        cmd = ['mafft']
        if self.threads > 1:
            cmd.extend(['--thread', str(self.threads)])
        cmd.append(input_path)
        
        with open(output_path, 'w') as output_file:
            result = subprocess.run(cmd, stdout=output_file, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"MAFFT failed: {result.stderr}")
    
    def _convert_clustal_to_fasta(self, clustal_path: str, fasta_path: str):
        """Convert ClustalW format to FASTA format."""
        sequences = {}
        
        with open(clustal_path, 'r') as f:
            lines = f.readlines()
        
        # Skip header lines
        start_idx = 0
        for i, line in enumerate(lines):
            if line.strip() == '' and i > 0:
                start_idx = i + 1
                break
        
        # Parse alignment blocks
        for line in lines[start_idx:]:
            line = line.strip()
            if line and not line.startswith('*') and not line.startswith(' '):
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    sequence = parts[1]
                    if seq_id not in sequences:
                        sequences[seq_id] = []
                    sequences[seq_id].append(sequence)
        
        # Write FASTA
        with open(fasta_path, 'w') as f:
            for seq_id, seq_parts in sequences.items():
                f.write(f">{seq_id}\n")
                full_sequence = ''.join(seq_parts)
                for i in range(0, len(full_sequence), 80):
                    f.write(f"{full_sequence[i:i+80]}\n")
    
    def _read_alignment_basic(self, file_path: str) -> List[SeqRecord]:
        """Basic alignment reader when BioPython is not available."""
        sequences = []
        
        with open(file_path, 'r') as f:
            current_id = None
            current_seq = []
            
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        seq_record = SeqRecord(
                            Seq(''.join(current_seq)),
                            id=current_id,
                            description=''
                        )
                        sequences.append(seq_record)
                    
                    current_id = line[1:].split()[0]
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
    
    def _basic_alignment(self, sequences: Dict[str, SeqRecord]) -> List[SeqRecord]:
        """
        Basic alignment fallback using simple gap insertion.
        This is a very basic implementation for when external tools are not available.
        """
        self.logger.warning("Using basic alignment - results may not be optimal")
        
        seq_list = list(sequences.values())
        
        # Find the longest sequence
        max_length = max(len(seq.seq) for seq in seq_list)
        
        # Simple gap padding to make all sequences the same length
        aligned_sequences = []
        for seq_record in seq_list:
            sequence = str(seq_record.seq)
            # Pad with gaps at the end
            padded_sequence = sequence + '-' * (max_length - len(sequence))
            
            aligned_record = SeqRecord(
                Seq(padded_sequence),
                id=seq_record.id,
                description=seq_record.description + " (basic alignment)"
            )
            aligned_sequences.append(aligned_record)
        
        return aligned_sequences
    
    def _save_alignment(self, alignment, output_path: str):
        """Save alignment to file."""
        if BIOPYTHON_AVAILABLE and hasattr(alignment, '_records'):
            AlignIO.write(alignment, output_path, 'fasta')
        else:
            # Basic save for list of SeqRecords
            with open(output_path, 'w') as f:
                if hasattr(alignment, '__iter__'):
                    sequence_list = alignment
                else:
                    sequence_list = [alignment]
                
                for seq_record in sequence_list:
                    f.write(f">{seq_record.id}\n")
                    sequence = str(seq_record.seq)
                    for i in range(0, len(sequence), 80):
                        f.write(f"{sequence[i:i+80]}\n")
        
        self.logger.info(f"Alignment saved to: {output_path}")
    
    def calculate_alignment_stats(self, alignment) -> Dict[str, float]:
        """
        Calculate alignment statistics.
        
        Args:
            alignment: MultipleSeqAlignment or list of SeqRecords
            
        Returns:
            Dictionary of alignment statistics
        """
        if hasattr(alignment, '_records'):
            sequences = [str(record.seq) for record in alignment]
        else:
            sequences = [str(record.seq) for record in alignment]
        
        if not sequences:
            return {}
        
        alignment_length = len(sequences[0])
        num_sequences = len(sequences)
        
        # Calculate gap percentage
        total_gaps = sum(seq.count('-') for seq in sequences)
        total_positions = alignment_length * num_sequences
        gap_percentage = (total_gaps / total_positions) * 100 if total_positions > 0 else 0
        
        # Calculate conservation
        conserved_positions = 0
        for pos in range(alignment_length):
            column = [seq[pos] for seq in sequences if pos < len(seq)]
            if len(set(column)) == 1 and column[0] != '-':
                conserved_positions += 1
        
        conservation_percentage = (conserved_positions / alignment_length) * 100 if alignment_length > 0 else 0
        
        stats = {
            'alignment_length': alignment_length,
            'num_sequences': num_sequences,
            'gap_percentage': gap_percentage,
            'conservation_percentage': conservation_percentage
        }
        
        return stats
    
    def trim_alignment(self, alignment, min_conservation: float = 0.5) -> List[SeqRecord]:
        """
        Trim alignment to remove poorly conserved regions.
        
        Args:
            alignment: MultipleSeqAlignment or list of SeqRecords
            min_conservation: Minimum conservation threshold (0-1)
            
        Returns:
            List of trimmed SeqRecord objects
        """
        if hasattr(alignment, '_records'):
            sequences = [str(record.seq) for record in alignment]
            records = list(alignment)
        else:
            sequences = [str(record.seq) for record in alignment]
            records = alignment
        
        if not sequences:
            return []
        
        alignment_length = len(sequences[0])
        num_sequences = len(sequences)
        
        # Identify positions to keep
        positions_to_keep = []
        
        for pos in range(alignment_length):
            column = [seq[pos] for seq in sequences if pos < len(seq)]
            non_gap_chars = [char for char in column if char != '-']
            
            # Calculate conservation at this position
            if non_gap_chars:
                most_common_char = max(set(non_gap_chars), key=non_gap_chars.count)
                conservation = non_gap_chars.count(most_common_char) / len(non_gap_chars)
                
                if conservation >= min_conservation:
                    positions_to_keep.append(pos)
        
        # Create trimmed sequences
        trimmed_records = []
        for i, record in enumerate(records):
            sequence = str(record.seq)
            trimmed_sequence = ''.join(sequence[pos] for pos in positions_to_keep if pos < len(sequence))
            
            trimmed_record = SeqRecord(
                Seq(trimmed_sequence),
                id=record.id,
                description=record.description + f" (trimmed, {len(positions_to_keep)}/{alignment_length} positions)"
            )
            trimmed_records.append(trimmed_record)
        
        self.logger.info(f"Trimmed alignment: {len(positions_to_keep)}/{alignment_length} positions retained")
        return trimmed_records