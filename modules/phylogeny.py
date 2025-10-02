"""
Phylogeny Module
===============

Handles phylogenetic tree construction and analysis.
Supports distance-based (NJ, UPGMA) and maximum likelihood methods.
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple, Any
from collections import defaultdict

try:
    from Bio import Phylo
    from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
    from Bio.Phylo import draw
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
try:
    import dendropy
    DENDROPY_AVAILABLE = True
except ImportError:
    DENDROPY_AVAILABLE = False

class PhylogeneticAnalyzer:
    """
    Phylogenetic tree construction and analysis.
    """
    
    def __init__(self, method: str = 'nj'):
        """
        Initialize phylogenetic analyzer.
        
        Args:
            method: Tree construction method ('nj', 'upgma', 'ml')
        """
        self.method = method.lower()
        self.logger = logging.getLogger(__name__)
        
        # Validate method
        if self.method not in ['nj', 'upgma', 'ml']:
            self.logger.warning(f"Unknown method '{method}'. Using 'nj' instead.")
            self.method = 'nj'
    
    def build_tree(self, alignment, output_path: Optional[str] = None):
        """
        Build phylogenetic tree from alignment.
        
        Args:
            alignment: MultipleSeqAlignment or list of SeqRecords
            output_path: Output file path for tree (Newick format)
            
        Returns:
            Tree object or Newick string
        """
        self.logger.info(f"Building phylogenetic tree using {self.method} method")
        
        if self.method == 'ml':
            return self._build_ml_tree(alignment, output_path)
        else:
            return self._build_distance_tree(alignment, output_path)
    
    def _build_distance_tree(self, alignment, output_path: Optional[str] = None):
        """Build distance-based tree (NJ or UPGMA)."""
        # Convert alignment to proper format
        if hasattr(alignment, '_records'):
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        else:
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        
        # Calculate distance matrix
        distance_matrix = self._calculate_distance_matrix(sequences, seq_ids)
        
        if BIOPYTHON_AVAILABLE:
            return self._build_biopython_tree(distance_matrix, output_path)
        else:
            return self._build_basic_tree(distance_matrix, seq_ids, output_path)
    
    def _calculate_distance_matrix(self, sequences: List[str], seq_ids: List[str]) -> np.ndarray:
        """Calculate pairwise distance matrix."""
        n = len(sequences)
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                distance = self._calculate_pairwise_distance(sequences[i], sequences[j])
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance
        
        return distance_matrix
    
    def _calculate_pairwise_distance(self, seq1: str, seq2: str) -> float:
        """Calculate pairwise distance between two sequences."""
        if len(seq1) != len(seq2):
            # Pad shorter sequence
            max_len = max(len(seq1), len(seq2))
            seq1 = seq1.ljust(max_len, '-')
            seq2 = seq2.ljust(max_len, '-')
        
        # Count differences (ignoring gaps)
        differences = 0
        valid_positions = 0
        
        for i in range(len(seq1)):
            if seq1[i] != '-' and seq2[i] != '-':
                valid_positions += 1
                if seq1[i] != seq2[i]:
                    differences += 1
        
        # Calculate p-distance (proportion of different sites)
        if valid_positions == 0:
            return 1.0  # Maximum distance
        
        p_distance = differences / valid_positions
        
        # Apply Jukes-Cantor correction for nucleotides
        # or simple p-distance for proteins
        if self._is_nucleotide_sequence(seq1):
            # Jukes-Cantor correction for nucleotides
            if p_distance < 0.75:
                return -0.75 * np.log(1 - (4/3) * p_distance)
            else:
                return 3.0  # Maximum corrected distance
        else:
            # Simple p-distance for proteins
            return p_distance
    
    def _is_nucleotide_sequence(self, sequence: str) -> bool:
        """Check if sequence is nucleotide or protein."""
        nucleotides = set('ATCGUN-')
        seq_chars = set(sequence.upper())
        return len(seq_chars - nucleotides) < len(seq_chars) * 0.1
    
    def _build_biopython_tree(self, distance_matrix: np.ndarray, output_path: Optional[str] = None):
        """Build tree using BioPython."""
        try:
            from Bio.Phylo.TreeConstruction import DistanceMatrix
            
            # Create sequence IDs
            seq_ids = [f"seq_{i}" for i in range(len(distance_matrix))]
            
            # Create BioPython DistanceMatrix object
            matrix_data = []
            for i in range(len(distance_matrix)):
                row = []
                for j in range(i):
                    row.append(distance_matrix[i][j])
                matrix_data.append(row)
            
            # Create DistanceMatrix object
            dm = DistanceMatrix(names=seq_ids, matrix=matrix_data)
            
            # Create tree constructor
            constructor = DistanceTreeConstructor()
            
            if self.method == 'nj':
                tree = constructor.nj(dm)
            elif self.method == 'upgma':
                tree = constructor.upgma(dm)
            
            # Save tree if output path provided
            if output_path:
                Phylo.write(tree, output_path, 'newick')
                self.logger.info(f"Tree saved to: {output_path}")
            
            return tree
            
        except Exception as e:
            self.logger.warning(f"BioPython tree construction failed: {e}. Using basic implementation.")
            return self._build_basic_tree(distance_matrix, [f"seq_{i}" for i in range(len(distance_matrix))], output_path)
    
    def _build_basic_tree(self, distance_matrix: np.ndarray, seq_ids: List[str], output_path: Optional[str] = None):
        """Build tree using basic implementation."""
        if self.method == 'nj':
            tree_newick = self._neighbor_joining(distance_matrix, seq_ids)
        elif self.method == 'upgma':
            tree_newick = self._upgma(distance_matrix, seq_ids)
        
        # Save tree if output path provided
        if output_path:
            with open(output_path, 'w') as f:
                f.write(tree_newick)
            self.logger.info(f"Tree saved to: {output_path}")
        
        return tree_newick
    
    def _neighbor_joining(self, distance_matrix: np.ndarray, seq_ids: List[str]) -> str:
        """Basic neighbor-joining implementation."""
        n = len(distance_matrix)
        if n < 2:
            return "();"
        
        # Copy matrix and sequence IDs
        D = distance_matrix.copy()
        taxa = seq_ids.copy()
        
        # Build tree iteratively
        while len(taxa) > 2:
            n = len(taxa)
            
            # Calculate Q matrix
            Q = np.zeros((n, n))
            for i in range(n):
                for j in range(n):
                    if i != j:
                        Q[i][j] = (n - 2) * D[i][j] - sum(D[i]) - sum(D[j])
            
            # Find minimum Q value
            min_val = float('inf')
            min_i, min_j = 0, 1
            for i in range(n):
                for j in range(i + 1, n):
                    if Q[i][j] < min_val:
                        min_val = Q[i][j]
                        min_i, min_j = i, j
            
            # Calculate branch lengths
            dist_i = (D[min_i][min_j] + (sum(D[min_i]) - sum(D[min_j])) / (n - 2)) / 2
            dist_j = D[min_i][min_j] - dist_i
            
            # Create new node
            new_taxon = f"({taxa[min_i]}:{dist_i:.6f},{taxa[min_j]}:{dist_j:.6f})"
            
            # Update distance matrix
            new_D = np.zeros((n - 1, n - 1))
            new_taxa = []
            
            # Add distances to new internal node
            k = 0
            for idx in range(n):
                if idx != min_i and idx != min_j:
                    dist_to_new = (D[min_i][idx] + D[min_j][idx] - D[min_i][min_j]) / 2
                    new_D[0][k + 1] = dist_to_new
                    new_D[k + 1][0] = dist_to_new
                    new_taxa.append(taxa[idx])
                    k += 1
            
            # Add other distances
            k1 = 1
            for i in range(n):
                if i != min_i and i != min_j:
                    k2 = 1
                    for j in range(n):
                        if j != min_i and j != min_j and j != i:
                            new_D[k1][k2] = D[i][j]
                            k2 += 1
                    k1 += 1
            
            # Update for next iteration
            D = new_D
            taxa = [new_taxon] + new_taxa
        
        # Final join
        if len(taxa) == 2:
            final_dist = D[0][1] / 2
            return f"({taxa[0]}:{final_dist:.6f},{taxa[1]}:{final_dist:.6f});"
        
        return f"({taxa[0]});"
    
    def _upgma(self, distance_matrix: np.ndarray, seq_ids: List[str]) -> str:
        """Basic UPGMA implementation."""
        n = len(distance_matrix)
        if n < 2:
            return "();"
        
        # Initialize clusters
        clusters = [[i] for i in range(n)]
        cluster_names = seq_ids.copy()
        D = distance_matrix.copy()
        heights = [0.0] * n
        
        while len(clusters) > 1:
            # Find minimum distance
            min_dist = float('inf')
            min_i, min_j = 0, 1
            
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    if D[i][j] < min_dist:
                        min_dist = D[i][j]
                        min_i, min_j = i, j
            
            # Calculate height of new node
            new_height = min_dist / 2
            
            # Calculate branch lengths
            branch_i = new_height - heights[min_i] if min_i < len(heights) else new_height
            branch_j = new_height - heights[min_j] if min_j < len(heights) else new_height
            
            # Create new cluster name
            new_name = f"({cluster_names[min_i]}:{branch_i:.6f},{cluster_names[min_j]}:{branch_j:.6f})"
            
            # Merge clusters
            new_cluster = clusters[min_i] + clusters[min_j]
            
            # Update distance matrix using UPGMA rule
            new_size = len(clusters) - 1
            new_D = np.zeros((new_size, new_size))
            new_clusters = []
            new_names = []
            new_heights_list = []
            
            # Add merged cluster
            new_clusters.append(new_cluster)
            new_names.append(new_name)
            new_heights_list.append(new_height)
            
            # Add other clusters
            idx = 1
            for k in range(len(clusters)):
                if k != min_i and k != min_j:
                    new_clusters.append(clusters[k])
                    new_names.append(cluster_names[k])
                    new_heights_list.append(heights[k] if k < len(heights) else 0.0)
                    
                    # Calculate distance to merged cluster
                    size_i = len(clusters[min_i])
                    size_j = len(clusters[min_j])
                    dist_merged = (size_i * D[min_i][k] + size_j * D[min_j][k]) / (size_i + size_j)
                    
                    new_D[0][idx] = dist_merged
                    new_D[idx][0] = dist_merged
                    idx += 1
            
            # Fill in other distances
            idx1 = 1
            for i in range(len(clusters)):
                if i != min_i and i != min_j:
                    idx2 = 1
                    for j in range(len(clusters)):
                        if j != min_i and j != min_j and j != i:
                            new_D[idx1][idx2] = D[i][j]
                            idx2 += 1
                    idx1 += 1
            
            # Update for next iteration
            clusters = new_clusters
            cluster_names = new_names
            heights = new_heights_list
            D = new_D
        
        return f"{cluster_names[0]};"
    
    def _build_ml_tree(self, alignment, output_path: Optional[str] = None):
        """Build maximum likelihood tree (basic implementation)."""
        self.logger.warning("Maximum likelihood method not fully implemented. Using NJ instead.")
        # Fall back to neighbor-joining for now
        old_method = self.method
        self.method = 'nj'
        result = self._build_distance_tree(alignment, output_path)
        self.method = old_method
        return result
    
    def visualize_tree(self, tree, output_path: Optional[str] = None, title: str = "Phylogenetic Tree"):
        """
        Visualize phylogenetic tree.
        
        Args:
            tree: Tree object or Newick string
            output_path: Output file path for plot
            title: Plot title
        """
        plt.figure(figsize=(12, 8))
        
        if BIOPYTHON_AVAILABLE and hasattr(tree, 'clade'):
            # BioPython tree
            Phylo.draw(tree, do_show=False)
        else:
            # Basic tree visualization for Newick string
            self._draw_basic_tree(tree)
        
        plt.title(title)
        plt.tight_layout()
        
        if output_path:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Tree plot saved to: {output_path}")
        else:
            plt.show()
        
        plt.close()
    
    def _draw_basic_tree(self, newick_string: str):
        """Basic tree drawing for Newick strings."""
        # This is a very simplified tree drawing
        # In practice, you would want to use a proper tree visualization library
        
        plt.text(0.5, 0.5, "Tree visualization requires BioPython or specialized libraries", 
                ha='center', va='center', fontsize=12, 
                bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue"))
        
        # Add the Newick string as text
        plt.text(0.5, 0.3, f"Newick: {newick_string[:100]}{'...' if len(newick_string) > 100 else ''}", 
                ha='center', va='center', fontsize=8, family='monospace')
        
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.axis('off')
    
    def calculate_tree_stats(self, tree) -> Dict[str, Any]:
        """
        Calculate tree statistics.
        
        Args:
            tree: Tree object or Newick string
            
        Returns:
            Dictionary of tree statistics
        """
        stats = {}
        
        if isinstance(tree, str):
            # Newick string - basic stats
            stats['tree_length'] = len(tree)
            stats['num_taxa'] = tree.count(',') + 1 if ',' in tree else 1
            stats['num_internal_nodes'] = tree.count('(')
        elif BIOPYTHON_AVAILABLE and hasattr(tree, 'clade'):
            # BioPython tree
            stats['num_taxa'] = tree.count_terminals()
            stats['tree_depth'] = tree.depths().get(tree.clade, 0)
            stats['total_branch_length'] = tree.total_branch_length()
        else:
            stats['note'] = "Tree statistics require proper tree object"
        
        return stats
    
    def get_clusters_from_tree(self, tree, threshold: float = 0.5) -> Dict[str, List[str]]:
        """
        Extract clusters from tree based on distance threshold.
        
        Args:
            tree: Tree object or Newick string
            threshold: Distance threshold for clustering
            
        Returns:
            Dictionary of cluster ID to list of sequence IDs
        """
        clusters = defaultdict(list)
        
        if isinstance(tree, str):
            # Basic clustering from Newick string
            # This is a simplified implementation
            cluster_id = 0
            # Extract terminal node names
            import re
            terminals = re.findall(r'([^(),:\s]+):', tree)
            
            # For now, put all terminals in one cluster
            # A proper implementation would parse the tree structure
            clusters[f"cluster_{cluster_id}"] = terminals
        
        elif BIOPYTHON_AVAILABLE and hasattr(tree, 'clade'):
            # BioPython tree clustering
            cluster_id = 0
            
            def traverse_and_cluster(clade, current_distance=0):
                nonlocal cluster_id
                
                if clade.is_terminal():
                    clusters[f"cluster_{cluster_id}"].append(clade.name)
                else:
                    if current_distance >= threshold:
                        # Start new cluster
                        cluster_id += 1
                    
                    for child in clade.clades:
                        traverse_and_cluster(child, current_distance + (child.branch_length or 0))
            
            traverse_and_cluster(tree.clade)
        
        return dict(clusters)