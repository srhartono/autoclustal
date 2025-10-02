"""
Clustering Module
================

Handles sequence clustering and representative selection.
Supports distance-based clustering and phylogenetic tree-based clustering.
"""

import logging
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from collections import defaultdict

try:
    from sklearn.cluster import AgglomerativeClustering, DBSCAN, KMeans
    from sklearn.metrics import silhouette_score
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

try:
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
class SequenceClusterer:
    """
    Sequence clustering and representative selection.
    """
    
    def __init__(self, threshold: float = 0.8, method: str = 'distance'):
        """
        Initialize sequence clusterer.
        
        Args:
            threshold: Clustering threshold (distance or similarity)
            method: Clustering method ('distance', 'tree', 'kmeans')
        """
        self.threshold = threshold
        self.method = method.lower()
        self.logger = logging.getLogger(__name__)
        
        if not SKLEARN_AVAILABLE and self.method in ['kmeans', 'dbscan']:
            self.logger.warning(f"scikit-learn not available. Using distance-based clustering.")
            self.method = 'distance'
    
    def cluster_sequences(self, alignment, tree=None) -> Dict[str, List[str]]:
        """
        Cluster sequences based on alignment and/or phylogenetic tree.
        
        Args:
            alignment: MultipleSeqAlignment or list of SeqRecords
            tree: Phylogenetic tree (optional)
            
        Returns:
            Dictionary of cluster ID to list of sequence IDs
        """
        self.logger.info(f"Clustering sequences using {self.method} method")
        
        # Extract sequences and IDs
        if hasattr(alignment, '_records'):
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        else:
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        
        if len(sequences) < 2:
            return {'cluster_0': seq_ids}
        
        # Choose clustering method
        if self.method == 'tree' and tree is not None:
            return self._cluster_by_tree(tree, seq_ids)
        elif self.method == 'kmeans':
            return self._cluster_by_kmeans(sequences, seq_ids)
        elif self.method == 'dbscan':
            return self._cluster_by_dbscan(sequences, seq_ids)
        else:
            return self._cluster_by_distance(sequences, seq_ids)
    
    def _calculate_distance_matrix(self, sequences: List[str]) -> np.ndarray:
        """Calculate pairwise distance matrix between sequences."""
        n = len(sequences)
        distance_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                distance = self._calculate_sequence_distance(sequences[i], sequences[j])
                distance_matrix[i][j] = distance
                distance_matrix[j][i] = distance
        
        return distance_matrix
    
    def _calculate_sequence_distance(self, seq1: str, seq2: str) -> float:
        """Calculate distance between two sequences."""
        if len(seq1) != len(seq2):
            # Align sequences to same length
            max_len = max(len(seq1), len(seq2))
            seq1 = seq1.ljust(max_len, '-')
            seq2 = seq2.ljust(max_len, '-')
        
        # Count mismatches (ignoring gaps)
        mismatches = 0
        valid_positions = 0
        
        for i in range(len(seq1)):
            if seq1[i] != '-' and seq2[i] != '-':
                valid_positions += 1
                if seq1[i] != seq2[i]:
                    mismatches += 1
        
        if valid_positions == 0:
            return 1.0  # Maximum distance
        
        return mismatches / valid_positions
    
    def _cluster_by_distance(self, sequences: List[str], seq_ids: List[str]) -> Dict[str, List[str]]:
        """Cluster sequences based on distance matrix."""
        distance_matrix = self._calculate_distance_matrix(sequences)
        
        if SKLEARN_AVAILABLE:
            return self._sklearn_hierarchical_clustering(distance_matrix, seq_ids)
        else:
            return self._basic_hierarchical_clustering(distance_matrix, seq_ids)
    
    def _sklearn_hierarchical_clustering(self, distance_matrix: np.ndarray, seq_ids: List[str]) -> Dict[str, List[str]]:
        """Use scikit-learn for hierarchical clustering."""
        # Convert distance to similarity for threshold
        similarity_threshold = 1.0 - self.threshold
        
        try:
            clustering = AgglomerativeClustering(
                n_clusters=None,
                distance_threshold=similarity_threshold,
                metric='precomputed',
                linkage='average'
            )
        except TypeError:
            # Fallback for older scikit-learn versions
            clustering = AgglomerativeClustering(
                n_clusters=None,
                distance_threshold=similarity_threshold,
                linkage='average'
            )
            # Use basic clustering if precomputed not supported
            return self._basic_hierarchical_clustering(distance_matrix, seq_ids)
        
        cluster_labels = clustering.fit_predict(distance_matrix)
        
        # Group sequences by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            clusters[f"cluster_{label}"].append(seq_ids[i])
        
        return dict(clusters)
    
    def _basic_hierarchical_clustering(self, distance_matrix: np.ndarray, seq_ids: List[str]) -> Dict[str, List[str]]:
        """Basic hierarchical clustering implementation."""
        n = len(seq_ids)
        
        # Initialize each sequence as its own cluster
        clusters = {i: [seq_ids[i]] for i in range(n)}
        cluster_distances = distance_matrix.copy()
        
        cluster_id = n
        
        while len(clusters) > 1:
            # Find minimum distance between clusters
            min_dist = float('inf')
            merge_i, merge_j = -1, -1
            
            cluster_indices = list(clusters.keys())
            for i in range(len(cluster_indices)):
                for j in range(i + 1, len(cluster_indices)):
                    idx_i, idx_j = cluster_indices[i], cluster_indices[j]
                    if idx_i < len(cluster_distances) and idx_j < len(cluster_distances):
                        dist = cluster_distances[idx_i][idx_j]
                        if dist < min_dist:
                            min_dist = dist
                            merge_i, merge_j = idx_i, idx_j
            
            # Stop if minimum distance exceeds threshold
            if min_dist > self.threshold:
                break
            
            # Merge clusters
            if merge_i in clusters and merge_j in clusters:
                new_cluster = clusters[merge_i] + clusters[merge_j]
                clusters[cluster_id] = new_cluster
                
                # Remove merged clusters
                del clusters[merge_i]
                del clusters[merge_j]
                
                cluster_id += 1
        
        # Rename clusters
        final_clusters = {}
        for i, (_, members) in enumerate(clusters.items()):
            final_clusters[f"cluster_{i}"] = members
        
        return final_clusters
    
    def _cluster_by_kmeans(self, sequences: List[str], seq_ids: List[str]) -> Dict[str, List[str]]:
        """Cluster sequences using K-means."""
        if not SKLEARN_AVAILABLE:
            return self._cluster_by_distance(sequences, seq_ids)
        
        # Convert sequences to feature vectors
        feature_matrix = self._sequences_to_features(sequences)
        
        # Determine optimal number of clusters
        max_clusters = min(10, len(sequences) // 2)
        if max_clusters < 2:
            return {'cluster_0': seq_ids}
        
        best_k = 2
        best_score = -1
        
        for k in range(2, max_clusters + 1):
            kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
            cluster_labels = kmeans.fit_predict(feature_matrix)
            
            if len(set(cluster_labels)) > 1:  # More than one cluster
                score = silhouette_score(feature_matrix, cluster_labels)
                if score > best_score:
                    best_score = score
                    best_k = k
        
        # Final clustering with best k
        kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(feature_matrix)
        
        # Group sequences by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            clusters[f"cluster_{label}"].append(seq_ids[i])
        
        return dict(clusters)
    
    def _cluster_by_dbscan(self, sequences: List[str], seq_ids: List[str]) -> Dict[str, List[str]]:
        """Cluster sequences using DBSCAN."""
        if not SKLEARN_AVAILABLE:
            return self._cluster_by_distance(sequences, seq_ids)
        
        # Convert sequences to feature vectors
        feature_matrix = self._sequences_to_features(sequences)
        
        # Use DBSCAN with distance threshold
        eps = 1.0 - self.threshold  # Convert similarity to distance
        dbscan = DBSCAN(eps=eps, min_samples=2, metric='euclidean')
        cluster_labels = dbscan.fit_predict(feature_matrix)
        
        # Group sequences by cluster
        clusters = defaultdict(list)
        for i, label in enumerate(cluster_labels):
            if label == -1:  # Noise points
                clusters[f"singleton_{i}"].append(seq_ids[i])
            else:
                clusters[f"cluster_{label}"].append(seq_ids[i])
        
        return dict(clusters)
    
    def _sequences_to_features(self, sequences: List[str]) -> np.ndarray:
        """Convert sequences to feature vectors for machine learning."""
        # Simple k-mer frequency features
        k = 3  # Tri-nucleotide/amino acid
        
        # Get all possible k-mers
        if self._is_nucleotide_sequence(sequences[0]):
            alphabet = 'ATCG'
        else:
            alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        
        kmers = []
        for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                for l in range(len(alphabet)):
                    kmers.append(alphabet[i] + alphabet[j] + alphabet[l])
        
        # Calculate k-mer frequencies for each sequence
        feature_matrix = np.zeros((len(sequences), len(kmers)))
        
        for seq_idx, sequence in enumerate(sequences):
            # Remove gaps
            clean_seq = sequence.replace('-', '')
            
            # Count k-mers
            kmer_counts = {}
            for i in range(len(clean_seq) - k + 1):
                kmer = clean_seq[i:i+k]
                if all(c in alphabet for c in kmer):
                    kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Normalize by sequence length
            total_kmers = sum(kmer_counts.values())
            if total_kmers > 0:
                for kmer_idx, kmer in enumerate(kmers):
                    feature_matrix[seq_idx][kmer_idx] = kmer_counts.get(kmer, 0) / total_kmers
        
        return feature_matrix
    
    def _is_nucleotide_sequence(self, sequence: str) -> bool:
        """Check if sequence is nucleotide or protein."""
        nucleotides = set('ATCGUN-')
        seq_chars = set(sequence.upper())
        return len(seq_chars - nucleotides) < len(seq_chars) * 0.1
    
    def _cluster_by_tree(self, tree, seq_ids: List[str]) -> Dict[str, List[str]]:
        """Cluster sequences based on phylogenetic tree."""
        # This is a simplified implementation
        # A more sophisticated approach would traverse the tree and cut at specific branch lengths
        
        if isinstance(tree, str):
            # Basic clustering from Newick string
            return self._cluster_from_newick(tree, seq_ids)
        elif hasattr(tree, 'clade'):
            # BioPython tree
            return self._cluster_from_biopython_tree(tree, seq_ids)
        else:
            # Fallback to distance-based clustering
            self.logger.warning("Unknown tree format. Using distance-based clustering.")
            return {'cluster_0': seq_ids}
    
    def _cluster_from_newick(self, newick_string: str, seq_ids: List[str]) -> Dict[str, List[str]]:
        """Extract clusters from Newick string (simplified)."""
        # This is a very basic implementation
        # For now, put all sequences in one cluster
        return {'cluster_0': seq_ids}
    
    def _cluster_from_biopython_tree(self, tree, seq_ids: List[str]) -> Dict[str, List[str]]:
        """Extract clusters from BioPython tree."""
        clusters = defaultdict(list)
        cluster_id = 0
        
        def traverse_and_cluster(clade, current_distance=0):
            nonlocal cluster_id
            
            if clade.is_terminal():
                clusters[f"cluster_{cluster_id}"].append(clade.name)
            else:
                # Check if we should split here based on branch length
                if current_distance >= self.threshold:
                    cluster_id += 1
                
                for child in clade.clades:
                    traverse_and_cluster(child, current_distance + (child.branch_length or 0))
        
        traverse_and_cluster(tree.clade)
        return dict(clusters)
    
    def select_representatives(self, clusters: Dict[str, List[str]], 
                            sequences: Dict[str, Any]) -> Dict[str, Any]:
        """
        Select representative sequences from each cluster.
        
        Args:
            clusters: Dictionary of cluster ID to list of sequence IDs
            sequences: Dictionary of sequence ID to SeqRecord objects
            
        Returns:
            Dictionary of cluster ID to representative SeqRecord
        """
        self.logger.info("Selecting representative sequences from clusters")
        
        representatives = {}
        
        for cluster_id, seq_ids in clusters.items():
            if len(seq_ids) == 1:
                # Single sequence cluster
                representatives[cluster_id] = sequences[seq_ids[0]]
            else:
                # Multiple sequences - select representative
                rep_seq = self._select_cluster_representative(seq_ids, sequences)
                representatives[cluster_id] = rep_seq
        
        self.logger.info(f"Selected {len(representatives)} representative sequences")
        return representatives
    
    def _select_cluster_representative(self, seq_ids: List[str], 
                                    sequences: Dict[str, Any]) -> Any:
        """
        Select the best representative from a cluster.
        
        Args:
            seq_ids: List of sequence IDs in the cluster
            sequences: Dictionary of sequence ID to SeqRecord objects
            
        Returns:
            Representative SeqRecord
        """
        if len(seq_ids) == 1:
            return sequences[seq_ids[0]]
        
        # Strategy: select the sequence closest to the centroid
        cluster_sequences = [str(sequences[seq_id].seq) for seq_id in seq_ids]
        
        # Calculate pairwise distances within cluster
        distances = []
        for i, seq1 in enumerate(cluster_sequences):
            total_distance = 0
            for j, seq2 in enumerate(cluster_sequences):
                if i != j:
                    total_distance += self._calculate_sequence_distance(seq1, seq2)
            distances.append(total_distance / (len(cluster_sequences) - 1))
        
        # Select sequence with minimum average distance (most central)
        min_dist_idx = distances.index(min(distances))
        representative_id = seq_ids[min_dist_idx]
        
        return sequences[representative_id]
    
    def calculate_cluster_stats(self, clusters: Dict[str, List[str]], 
                              sequences: Dict[str, Any]) -> Dict[str, Dict[str, float]]:
        """
        Calculate statistics for each cluster.
        
        Args:
            clusters: Dictionary of cluster ID to list of sequence IDs
            sequences: Dictionary of sequence ID to SeqRecord objects
            
        Returns:
            Dictionary of cluster statistics
        """
        stats = {}
        
        for cluster_id, seq_ids in clusters.items():
            cluster_sequences = [str(sequences[seq_id].seq) for seq_id in seq_ids]
            
            if len(cluster_sequences) == 1:
                cluster_stats = {
                    'size': 1,
                    'avg_length': len(cluster_sequences[0]),
                    'intra_cluster_distance': 0.0,
                    'conservation': 1.0
                }
            else:
                # Calculate intra-cluster distances
                distances = []
                for i in range(len(cluster_sequences)):
                    for j in range(i + 1, len(cluster_sequences)):
                        dist = self._calculate_sequence_distance(
                            cluster_sequences[i], cluster_sequences[j]
                        )
                        distances.append(dist)
                
                # Calculate conservation (1 - average distance)
                avg_distance = sum(distances) / len(distances) if distances else 0
                conservation = 1.0 - avg_distance
                
                cluster_stats = {
                    'size': len(seq_ids),
                    'avg_length': sum(len(seq.replace('-', '')) for seq in cluster_sequences) / len(cluster_sequences),
                    'intra_cluster_distance': avg_distance,
                    'conservation': conservation
                }
            
            stats[cluster_id] = cluster_stats
        
        return stats
    
    def optimize_clustering_threshold(self, sequences: List[str], seq_ids: List[str]) -> float:
        """
        Optimize clustering threshold using silhouette analysis.
        
        Args:
            sequences: List of sequence strings
            seq_ids: List of sequence IDs
            
        Returns:
            Optimized threshold value
        """
        if not SKLEARN_AVAILABLE:
            self.logger.warning("scikit-learn not available. Using default threshold.")
            return self.threshold
        
        distance_matrix = self._calculate_distance_matrix(sequences)
        
        # Try different thresholds
        thresholds = np.arange(0.1, 1.0, 0.1)
        best_threshold = self.threshold
        best_score = -1
        
        for threshold in thresholds:
            old_threshold = self.threshold
            self.threshold = threshold
            
            try:
                clustering = AgglomerativeClustering(
                    n_clusters=None,
                    distance_threshold=threshold,
                    affinity='precomputed',
                    linkage='average'
                )
                
                cluster_labels = clustering.fit_predict(distance_matrix)
                
                # Check if we have more than one cluster
                if len(set(cluster_labels)) > 1:
                    score = silhouette_score(distance_matrix, cluster_labels, metric='precomputed')
                    if score > best_score:
                        best_score = score
                        best_threshold = threshold
            
            except Exception as e:
                self.logger.debug(f"Failed to evaluate threshold {threshold}: {e}")
            
            finally:
                self.threshold = old_threshold
        
        self.logger.info(f"Optimized clustering threshold: {best_threshold}")
        return best_threshold