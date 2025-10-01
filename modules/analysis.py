"""
Analysis and Reporting Module
============================

Handles comprehensive analysis and reporting of sequence analysis results.
Includes PCA analysis, summary tables, and visualization.
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.manifold import TSNE
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

class AnalysisReporter:
    """
    Comprehensive analysis and reporting for sequence analysis results.
    """
    
    def __init__(self, output_dir: Path):
        """
        Initialize analysis reporter.
        
        Args:
            output_dir: Output directory for reports and plots
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
        
        # Set up matplotlib style
        plt.style.use('default')
        sns.set_palette("husl")
    
    def generate_summary_table(self, sequences: Dict[str, Any], 
                             clusters: Dict[str, List[str]], 
                             annotations: Dict[str, Any] = None):
        """
        Generate comprehensive summary table.
        
        Args:
            sequences: Dictionary of sequence ID to SeqRecord
            clusters: Dictionary of cluster ID to list of sequence IDs
            annotations: Dictionary of sequence annotations from BLAST/BLAT
        """
        self.logger.info("Generating summary table")
        
        # Prepare data for summary table
        summary_data = []
        
        for cluster_id, seq_ids in clusters.items():
            for seq_id in seq_ids:
                if seq_id in sequences:
                    seq_record = sequences[seq_id]
                    
                    # Basic sequence information
                    row = {
                        'sequence_id': seq_id,
                        'cluster_id': cluster_id,
                        'sequence_length': len(str(seq_record.seq).replace('-', '')),
                        'description': getattr(seq_record, 'description', ''),
                    }
                    
                    # Add annotation information if available
                    if annotations and seq_id in annotations:
                        annotation = annotations[seq_id]
                        if 'hits' in annotation and annotation['hits']:
                            best_hit = annotation['hits'][0]
                            row.update({
                                'best_match_title': best_hit.get('title', ''),
                                'best_match_accession': best_hit.get('accession', ''),
                                'best_match_evalue': best_hit.get('evalue', ''),
                                'best_match_identity': best_hit.get('identity', ''),
                                'organism': self._extract_organism_from_hit(best_hit)
                            })
                        else:
                            row.update({
                                'best_match_title': 'No significant hits',
                                'best_match_accession': '',
                                'best_match_evalue': '',
                                'best_match_identity': '',
                                'organism': 'Unknown'
                            })
                    else:
                        row.update({
                            'best_match_title': 'Not searched',
                            'best_match_accession': '',
                            'best_match_evalue': '',
                            'best_match_identity': '',
                            'organism': 'Unknown'
                        })
                    
                    summary_data.append(row)
        
        # Create summary table
        if PANDAS_AVAILABLE:
            df = pd.DataFrame(summary_data)
            self._save_dataframe_reports(df, 'summary_table')
        else:
            self._save_basic_table(summary_data, 'summary_table.txt')
        
        # Generate cluster statistics
        self._generate_cluster_statistics(clusters, sequences, annotations)
    
    def _extract_organism_from_hit(self, hit: Dict[str, Any]) -> str:
        """Extract organism name from BLAST hit."""
        title = hit.get('title', '')
        
        # Try to extract organism from [brackets]
        import re
        bracket_match = re.search(r'\[([^\]]+)\]', title)
        if bracket_match:
            return bracket_match.group(1)
        
        # Try common patterns
        for pattern in ['OS=', 'organism=']:
            if pattern in title.lower():
                parts = title.lower().split(pattern)
                if len(parts) > 1:
                    org_words = parts[1].split()[:2]  # Take first two words
                    return ' '.join(org_words)
        
        return 'Unknown'
    
    def _save_dataframe_reports(self, df: 'pd.DataFrame', base_name: str):
        """Save DataFrame in multiple formats."""
        # Save as CSV
        csv_path = self.output_dir / f"{base_name}.csv"
        df.to_csv(csv_path, index=False)
        self.logger.info(f"Summary table saved as CSV: {csv_path}")
        
        # Save as Excel if possible
        try:
            excel_path = self.output_dir / f"{base_name}.xlsx"
            df.to_excel(excel_path, index=False)
            self.logger.info(f"Summary table saved as Excel: {excel_path}")
        except ImportError:
            self.logger.warning("openpyxl not available. Skipping Excel export.")
        
        # Save as HTML
        html_path = self.output_dir / f"{base_name}.html"
        html_content = df.to_html(index=False, escape=False)
        with open(html_path, 'w') as f:
            f.write(f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>AutoClustal Summary Table</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 20px; }}
                    table {{ border-collapse: collapse; width: 100%; }}
                    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                    th {{ background-color: #f2f2f2; font-weight: bold; }}
                    tr:nth-child(even) {{ background-color: #f9f9f9; }}
                </style>
            </head>
            <body>
                <h1>AutoClustal Analysis Summary</h1>
                {html_content}
            </body>
            </html>
            """)
        self.logger.info(f"Summary table saved as HTML: {html_path}")
    
    def _save_basic_table(self, data: List[Dict], filename: str):
        """Save table in basic text format when pandas is not available."""
        filepath = self.output_dir / filename
        
        if not data:
            with open(filepath, 'w') as f:
                f.write("No data to report\n")
            return
        
        # Get all keys
        all_keys = set()
        for row in data:
            all_keys.update(row.keys())
        
        headers = sorted(all_keys)
        
        with open(filepath, 'w') as f:
            # Write header
            f.write('\t'.join(headers) + '\n')
            
            # Write data rows
            for row in data:
                values = [str(row.get(key, '')) for key in headers]
                f.write('\t'.join(values) + '\n')
        
        self.logger.info(f"Summary table saved as text: {filepath}")
    
    def _generate_cluster_statistics(self, clusters: Dict[str, List[str]], 
                                   sequences: Dict[str, Any], 
                                   annotations: Dict[str, Any] = None):
        """Generate cluster-level statistics."""
        cluster_stats = []
        
        for cluster_id, seq_ids in clusters.items():
            # Calculate cluster size and sequence lengths
            cluster_sequences = [sequences[seq_id] for seq_id in seq_ids if seq_id in sequences]
            lengths = [len(str(seq.seq).replace('-', '')) for seq in cluster_sequences]
            
            # Collect organisms from annotations
            organisms = []
            if annotations:
                for seq_id in seq_ids:
                    if seq_id in annotations and 'hits' in annotations[seq_id]:
                        hits = annotations[seq_id]['hits']
                        if hits:
                            organism = self._extract_organism_from_hit(hits[0])
                            organisms.append(organism)
            
            # Most common organism
            if organisms:
                organism_counts = {}
                for org in organisms:
                    organism_counts[org] = organism_counts.get(org, 0) + 1
                most_common_organism = max(organism_counts.items(), key=lambda x: x[1])
                predominant_organism = f"{most_common_organism[0]} ({most_common_organism[1]}/{len(organisms)})"
            else:
                predominant_organism = "Unknown"
            
            stats = {
                'cluster_id': cluster_id,
                'cluster_size': len(seq_ids),
                'avg_sequence_length': np.mean(lengths) if lengths else 0,
                'min_sequence_length': min(lengths) if lengths else 0,
                'max_sequence_length': max(lengths) if lengths else 0,
                'predominant_organism': predominant_organism,
                'total_organisms': len(set(organisms)) if organisms else 0
            }
            
            cluster_stats.append(stats)
        
        # Save cluster statistics
        if PANDAS_AVAILABLE:
            cluster_df = pd.DataFrame(cluster_stats)
            self._save_dataframe_reports(cluster_df, 'cluster_statistics')
        else:
            self._save_basic_table(cluster_stats, 'cluster_statistics.txt')
    
    def perform_pca_analysis(self, alignment, clusters: Dict[str, List[str]], 
                           annotations: Dict[str, Any] = None):
        """
        Perform PCA analysis on sequence alignment.
        
        Args:
            alignment: Multiple sequence alignment
            clusters: Dictionary of clusters
            annotations: BLAST/BLAT annotations
        """
        self.logger.info("Performing PCA analysis")
        
        if not SKLEARN_AVAILABLE:
            self.logger.warning("scikit-learn not available. Skipping PCA analysis.")
            return
        
        # Extract sequences and metadata
        if hasattr(alignment, '_records'):
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        else:
            sequences = [str(record.seq) for record in alignment]
            seq_ids = [record.id for record in alignment]
        
        # Convert sequences to feature matrix
        feature_matrix = self._sequences_to_features(sequences)
        
        if feature_matrix.shape[1] < 2:
            self.logger.warning("Insufficient features for PCA analysis.")
            return
        
        # Standardize features
        scaler = StandardScaler()
        features_scaled = scaler.fit_transform(feature_matrix)
        
        # Perform PCA
        n_components = min(10, features_scaled.shape[1], features_scaled.shape[0])
        pca = PCA(n_components=n_components)
        pca_result = pca.fit_transform(features_scaled)
        
        # Prepare data for plotting
        pca_data = {
            'sequence_id': seq_ids,
            'PC1': pca_result[:, 0],
            'PC2': pca_result[:, 1] if pca_result.shape[1] > 1 else np.zeros(len(seq_ids))
        }
        
        # Add cluster information
        seq_to_cluster = {}
        for cluster_id, cluster_seq_ids in clusters.items():
            for seq_id in cluster_seq_ids:
                seq_to_cluster[seq_id] = cluster_id
        
        pca_data['cluster'] = [seq_to_cluster.get(seq_id, 'Unknown') for seq_id in seq_ids]
        
        # Add organism information if available
        if annotations:
            organisms = []
            for seq_id in seq_ids:
                if seq_id in annotations and 'hits' in annotations[seq_id]:
                    hits = annotations[seq_id]['hits']
                    if hits:
                        organism = self._extract_organism_from_hit(hits[0])
                        organisms.append(organism)
                    else:
                        organisms.append('Unknown')
                else:
                    organisms.append('Unknown')
            pca_data['organism'] = organisms
        else:
            pca_data['organism'] = ['Unknown'] * len(seq_ids)
        
        # Create PCA plots
        self._create_pca_plots(pca_data, pca.explained_variance_ratio_)
        
        # Perform t-SNE if enough samples
        if len(sequences) >= 10:
            self._perform_tsne_analysis(features_scaled, pca_data)
        
        # Save PCA data
        if PANDAS_AVAILABLE:
            pca_df = pd.DataFrame(pca_data)
            self._save_dataframe_reports(pca_df, 'pca_results')
    
    def _sequences_to_features(self, sequences: List[str]) -> np.ndarray:
        """Convert sequences to feature matrix for PCA."""
        # Use k-mer frequencies as features
        k = 3  # Use tri-nucleotides or tri-peptides
        
        # Determine if sequences are nucleotide or protein
        is_nucleotide = self._is_nucleotide_sequence(sequences[0])
        
        if is_nucleotide:
            alphabet = 'ATCG'
        else:
            alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        
        # Generate all possible k-mers
        kmers = []
        def generate_kmers(current, remaining):
            if remaining == 0:
                kmers.append(current)
                return
            for char in alphabet:
                generate_kmers(current + char, remaining - 1)
        
        generate_kmers('', k)
        
        # Calculate k-mer frequencies
        feature_matrix = np.zeros((len(sequences), len(kmers)))
        
        for seq_idx, sequence in enumerate(sequences):
            # Remove gaps
            clean_seq = sequence.replace('-', '').upper()
            
            # Count k-mers
            kmer_counts = {}
            for i in range(len(clean_seq) - k + 1):
                kmer = clean_seq[i:i+k]
                if all(c in alphabet for c in kmer):
                    kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Normalize by total k-mers
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
    
    def _create_pca_plots(self, pca_data: Dict[str, List], explained_variance: np.ndarray):
        """Create PCA visualization plots."""
        if PLOTLY_AVAILABLE:
            self._create_plotly_pca(pca_data, explained_variance)
        
        # Also create matplotlib plots
        self._create_matplotlib_pca(pca_data, explained_variance)
    
    def _create_plotly_pca(self, pca_data: Dict[str, List], explained_variance: np.ndarray):
        """Create interactive PCA plots using Plotly."""
        # Create subplot figure
        fig = make_subplots(
            rows=1, cols=2,
            subplot_titles=('PCA by Cluster', 'PCA by Organism'),
            specs=[[{'type': 'scatter'}, {'type': 'scatter'}]]
        )
        
        # Get unique clusters and organisms for colors
        unique_clusters = list(set(pca_data['cluster']))
        unique_organisms = list(set(pca_data['organism']))
        
        # Plot by cluster
        for i, cluster in enumerate(unique_clusters):
            mask = [c == cluster for c in pca_data['cluster']]
            x_vals = [pca_data['PC1'][j] for j, m in enumerate(mask) if m]
            y_vals = [pca_data['PC2'][j] for j, m in enumerate(mask) if m]
            names = [pca_data['sequence_id'][j] for j, m in enumerate(mask) if m]
            
            fig.add_trace(
                go.Scatter(
                    x=x_vals, y=y_vals,
                    mode='markers',
                    name=cluster,
                    text=names,
                    hovertemplate='<b>%{text}</b><br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>'
                ),
                row=1, col=1
            )
        
        # Plot by organism (limit to top 10 most common)
        organism_counts = {org: pca_data['organism'].count(org) for org in unique_organisms}
        top_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        for i, (organism, _) in enumerate(top_organisms):
            mask = [o == organism for o in pca_data['organism']]
            x_vals = [pca_data['PC1'][j] for j, m in enumerate(mask) if m]
            y_vals = [pca_data['PC2'][j] for j, m in enumerate(mask) if m]
            names = [pca_data['sequence_id'][j] for j, m in enumerate(mask) if m]
            
            fig.add_trace(
                go.Scatter(
                    x=x_vals, y=y_vals,
                    mode='markers',
                    name=organism,
                    text=names,
                    hovertemplate='<b>%{text}</b><br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>',
                    showlegend=False
                ),
                row=1, col=2
            )
        
        # Update layout
        fig.update_xaxes(title_text=f'PC1 ({explained_variance[0]:.1%} variance)', row=1, col=1)
        fig.update_xaxes(title_text=f'PC1 ({explained_variance[0]:.1%} variance)', row=1, col=2)
        fig.update_yaxes(title_text=f'PC2 ({explained_variance[1]:.1%} variance)', row=1, col=1)
        fig.update_yaxes(title_text=f'PC2 ({explained_variance[1]:.1%} variance)', row=1, col=2)
        
        fig.update_layout(
            title='PCA Analysis of Sequence Data',
            height=600,
            width=1200
        )
        
        # Save interactive plot
        html_path = self.output_dir / 'pca_analysis.html'
        fig.write_html(str(html_path))
        self.logger.info(f"Interactive PCA plot saved: {html_path}")
    
    def _create_matplotlib_pca(self, pca_data: Dict[str, List], explained_variance: np.ndarray):
        """Create static PCA plots using matplotlib."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot by cluster
        unique_clusters = list(set(pca_data['cluster']))
        colors = plt.cm.Set1(np.linspace(0, 1, len(unique_clusters)))
        
        for i, cluster in enumerate(unique_clusters):
            mask = [c == cluster for c in pca_data['cluster']]
            x_vals = [pca_data['PC1'][j] for j, m in enumerate(mask) if m]
            y_vals = [pca_data['PC2'][j] for j, m in enumerate(mask) if m]
            
            ax1.scatter(x_vals, y_vals, c=[colors[i]], label=cluster, alpha=0.7, s=50)
        
        ax1.set_xlabel(f'PC1 ({explained_variance[0]:.1%} variance)')
        ax1.set_ylabel(f'PC2 ({explained_variance[1]:.1%} variance)')
        ax1.set_title('PCA by Cluster')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # Plot by organism (top 10 most common)
        organism_counts = {org: pca_data['organism'].count(org) for org in set(pca_data['organism'])}
        top_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        colors_org = plt.cm.tab10(np.linspace(0, 1, len(top_organisms)))
        
        for i, (organism, _) in enumerate(top_organisms):
            mask = [o == organism for o in pca_data['organism']]
            x_vals = [pca_data['PC1'][j] for j, m in enumerate(mask) if m]
            y_vals = [pca_data['PC2'][j] for j, m in enumerate(mask) if m]
            
            ax2.scatter(x_vals, y_vals, c=[colors_org[i]], label=organism, alpha=0.7, s=50)
        
        ax2.set_xlabel(f'PC1 ({explained_variance[0]:.1%} variance)')
        ax2.set_ylabel(f'PC2 ({explained_variance[1]:.1%} variance)')
        ax2.set_title('PCA by Organism (Top 10)')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        png_path = self.output_dir / 'pca_analysis.png'
        plt.savefig(png_path, dpi=300, bbox_inches='tight')
        self.logger.info(f"PCA plot saved: {png_path}")
        
        plt.close()
    
    def _perform_tsne_analysis(self, features_scaled: np.ndarray, pca_data: Dict[str, List]):
        """Perform t-SNE analysis for non-linear dimensionality reduction."""
        try:
            # Perform t-SNE
            tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(features_scaled)-1))
            tsne_result = tsne.fit_transform(features_scaled)
            
            # Create t-SNE plot
            plt.figure(figsize=(12, 5))
            
            # Plot by cluster
            plt.subplot(1, 2, 1)
            unique_clusters = list(set(pca_data['cluster']))
            colors = plt.cm.Set1(np.linspace(0, 1, len(unique_clusters)))
            
            for i, cluster in enumerate(unique_clusters):
                mask = [c == cluster for c in pca_data['cluster']]
                x_vals = [tsne_result[j, 0] for j, m in enumerate(mask) if m]
                y_vals = [tsne_result[j, 1] for j, m in enumerate(mask) if m]
                
                plt.scatter(x_vals, y_vals, c=[colors[i]], label=cluster, alpha=0.7, s=50)
            
            plt.xlabel('t-SNE 1')
            plt.ylabel('t-SNE 2')
            plt.title('t-SNE by Cluster')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(True, alpha=0.3)
            
            # Plot by organism
            plt.subplot(1, 2, 2)
            organism_counts = {org: pca_data['organism'].count(org) for org in set(pca_data['organism'])}
            top_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:10]
            
            colors_org = plt.cm.tab10(np.linspace(0, 1, len(top_organisms)))
            
            for i, (organism, _) in enumerate(top_organisms):
                mask = [o == organism for o in pca_data['organism']]
                x_vals = [tsne_result[j, 0] for j, m in enumerate(mask) if m]
                y_vals = [tsne_result[j, 1] for j, m in enumerate(mask) if m]
                
                plt.scatter(x_vals, y_vals, c=[colors_org[i]], label=organism, alpha=0.7, s=50)
            
            plt.xlabel('t-SNE 1')
            plt.ylabel('t-SNE 2')
            plt.title('t-SNE by Organism (Top 10)')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            
            # Save t-SNE plot
            tsne_path = self.output_dir / 'tsne_analysis.png'
            plt.savefig(tsne_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"t-SNE plot saved: {tsne_path}")
            
            plt.close()
            
        except Exception as e:
            self.logger.warning(f"t-SNE analysis failed: {str(e)}")
    
    def create_visualizations(self, tree=None, clusters: Dict[str, List[str]] = None, 
                            annotations: Dict[str, Any] = None):
        """Create various visualizations of the analysis results."""
        self.logger.info("Creating visualizations")
        
        # Tree visualization
        if tree is not None:
            self._create_tree_visualization(tree)
        
        # Cluster size distribution
        if clusters:
            self._create_cluster_visualizations(clusters, annotations)
        
        # Annotation summary
        if annotations:
            self._create_annotation_visualizations(annotations)
    
    def _create_tree_visualization(self, tree):
        """Create phylogenetic tree visualization."""
        try:
            plt.figure(figsize=(12, 8))
            
            if hasattr(tree, 'clade'):  # BioPython tree
                from Bio import Phylo
                Phylo.draw(tree, do_show=False)
            else:
                # Basic tree representation for Newick strings
                plt.text(0.5, 0.5, "Phylogenetic Tree", ha='center', va='center', fontsize=16)
                plt.text(0.5, 0.4, "(Detailed tree visualization requires BioPython)", 
                        ha='center', va='center', fontsize=10)
                
                if isinstance(tree, str) and len(tree) < 500:
                    plt.text(0.5, 0.3, f"Newick: {tree}", ha='center', va='center', 
                            fontsize=8, family='monospace', wrap=True)
            
            plt.title('Phylogenetic Tree')
            
            # Save tree plot
            tree_path = self.output_dir / 'phylogenetic_tree.png'
            plt.savefig(tree_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Tree visualization saved: {tree_path}")
            
            plt.close()
            
        except Exception as e:
            self.logger.warning(f"Tree visualization failed: {str(e)}")
    
    def _create_cluster_visualizations(self, clusters: Dict[str, List[str]], 
                                     annotations: Dict[str, Any] = None):
        """Create cluster-related visualizations."""
        # Cluster size distribution
        cluster_sizes = [len(seq_ids) for seq_ids in clusters.values()]
        
        # Create cluster info with best match titles
        cluster_info = {}
        for cluster_id, seq_ids in clusters.items():
            cluster_info[cluster_id] = {
                'size': len(seq_ids),
                'best_match': 'No annotation'
            }
            
            # Try to get best match from annotations
            if annotations:
                best_matches = []
                for seq_id in seq_ids:
                    if seq_id in annotations and annotations[seq_id] and 'hits' in annotations[seq_id]:
                        hits = annotations[seq_id]['hits']
                        if hits:
                            title = hits[0].get('title', 'Unknown')
                            best_matches.append(title)
                
                if best_matches:
                    # Use the first best match (they should be similar within a cluster)
                    best_match = best_matches[0]
                    # Truncate long titles
                    if len(best_match) > 30:
                        best_match = best_match[:27] + "..."
                    cluster_info[cluster_id]['best_match'] = best_match
        
        plt.figure(figsize=(12, 8))
        
        plt.subplot(2, 1, 1)
        plt.hist(cluster_sizes, bins=max(1, len(set(cluster_sizes))), alpha=0.7, edgecolor='black')
        plt.xlabel('Cluster Size')
        plt.ylabel('Number of Clusters')
        plt.title('Cluster Size Distribution')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(2, 1, 2)
        cluster_names = list(clusters.keys())
        
        # Create labels with cluster ID and best match
        cluster_labels = []
        for cluster_id in cluster_names:
            best_match = cluster_info[cluster_id]['best_match']
            if best_match != 'No annotation':
                label = f"{cluster_id}\n{best_match}"
            else:
                label = cluster_id
            cluster_labels.append(label)
        
        bars = plt.barh(range(len(cluster_names)), cluster_sizes, color='skyblue', alpha=0.8, edgecolor='black')
        plt.yticks(range(len(cluster_names)), cluster_labels)
        plt.xlabel('Number of Sequences')
        plt.ylabel('Cluster ID + Best Match')
        plt.title('Sequences per Cluster (with Annotations)')
        plt.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for i, (bar, size) in enumerate(zip(bars, cluster_sizes)):
            plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                    str(size), ha='left', va='center', fontweight='bold')
        
        plt.tight_layout()
        
        # Save cluster plot
        cluster_path = self.output_dir / 'cluster_analysis.png'
        plt.savefig(cluster_path, dpi=300, bbox_inches='tight')
        self.logger.info(f"Cluster visualization saved: {cluster_path}")
        
        plt.close()
    
    def _create_annotation_visualizations(self, annotations: Dict[str, Any]):
        """Create annotation-related visualizations."""
        # Extract organism information
        organisms = []
        evalues = []
        identities = []
        
        for seq_id, annotation in annotations.items():
            if 'hits' in annotation and annotation['hits']:
                best_hit = annotation['hits'][0]
                organism = self._extract_organism_from_hit(best_hit)
                organisms.append(organism)
                
                if 'evalue' in best_hit:
                    evalues.append(best_hit['evalue'])
                
                if 'identity' in best_hit:
                    identities.append(best_hit['identity'])
        
        # Create organism distribution plot
        if organisms:
            organism_counts = {}
            for org in organisms:
                organism_counts[org] = organism_counts.get(org, 0) + 1
            
            # Top 15 most common organisms
            top_organisms = sorted(organism_counts.items(), key=lambda x: x[1], reverse=True)[:15]
            
            plt.figure(figsize=(12, 8))
            
            # Organism bar plot
            plt.subplot(2, 2, 1)
            orgs, counts = zip(*top_organisms)
            plt.barh(orgs, counts)
            plt.xlabel('Number of Sequences')
            plt.title('Top Organisms (BLAST Hits)')
            plt.gca().invert_yaxis()
            
            # E-value distribution
            if evalues:
                plt.subplot(2, 2, 2)
                plt.hist([np.log10(e) for e in evalues if e > 0], bins=20, alpha=0.7, edgecolor='black')
                plt.xlabel('log10(E-value)')
                plt.ylabel('Number of Hits')
                plt.title('E-value Distribution')
                plt.grid(True, alpha=0.3)
            
            # Identity distribution
            if identities:
                plt.subplot(2, 2, 3)
                plt.hist(identities, bins=20, alpha=0.7, edgecolor='black')
                plt.xlabel('Sequence Identity (%)')
                plt.ylabel('Number of Hits')
                plt.title('Identity Distribution')
                plt.grid(True, alpha=0.3)
            
            # Organism pie chart (top 10)
            if len(top_organisms) > 1:
                plt.subplot(2, 2, 4)
                top_10 = top_organisms[:10]
                others_count = sum(count for _, count in top_organisms[10:])
                
                pie_labels = [org for org, _ in top_10]
                pie_sizes = [count for _, count in top_10]
                
                if others_count > 0:
                    pie_labels.append('Others')
                    pie_sizes.append(others_count)
                
                plt.pie(pie_sizes, labels=pie_labels, autopct='%1.1f%%', startangle=90)
                plt.title('Organism Distribution')
            
            plt.tight_layout()
            
            # Save annotation plot
            annotation_path = self.output_dir / 'annotation_analysis.png'
            plt.savefig(annotation_path, dpi=300, bbox_inches='tight')
            self.logger.info(f"Annotation visualization saved: {annotation_path}")
            
            plt.close()
    
    def generate_final_report(self, sequences: Dict[str, Any], clusters: Dict[str, List[str]], 
                            annotations: Dict[str, Any] = None, tree=None):
        """Generate a comprehensive final report."""
        self.logger.info("Generating final report")
        
        report_path = self.output_dir / 'AutoClustal_Report.html'
        
        # Calculate summary statistics
        total_sequences = len(sequences)
        total_clusters = len(clusters)
        
        # Get cluster size statistics
        cluster_sizes = [len(seq_ids) for seq_ids in clusters.values()]
        avg_cluster_size = np.mean(cluster_sizes) if cluster_sizes else 0
        largest_cluster = max(cluster_sizes) if cluster_sizes else 0
        
        # Get sequence length statistics
        seq_lengths = [len(str(seq.seq).replace('-', '')) for seq in sequences.values()]
        avg_seq_length = np.mean(seq_lengths) if seq_lengths else 0
        
        # Get organism statistics if available
        organism_stats = ""
        if annotations:
            organisms = []
            for annotation in annotations.values():
                if 'hits' in annotation and annotation['hits']:
                    organism = self._extract_organism_from_hit(annotation['hits'][0])
                    organisms.append(organism)
            
            if organisms:
                unique_organisms = len(set(organisms))
                organism_counts = {}
                for org in organisms:
                    organism_counts[org] = organism_counts.get(org, 0) + 1
                
                top_organism = max(organism_counts.items(), key=lambda x: x[1])
                organism_stats = f"""
                <h3>Organism Analysis</h3>
                <ul>
                    <li>Total unique organisms identified: {unique_organisms}</li>
                    <li>Most common organism: {top_organism[0]} ({top_organism[1]} sequences)</li>
                    <li>Sequences with BLAST hits: {len(organisms)}/{total_sequences}</li>
                </ul>
                """
        
        # Generate HTML report
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>AutoClustal Analysis Report</title>
            <style>
                body {{
                    font-family: Arial, sans-serif;
                    line-height: 1.6;
                    margin: 40px;
                    background-color: #f5f5f5;
                }}
                .container {{
                    background-color: white;
                    padding: 30px;
                    border-radius: 10px;
                    box-shadow: 0 2px 10px rgba(0,0,0,0.1);
                }}
                h1 {{
                    color: #2c3e50;
                    border-bottom: 3px solid #3498db;
                    padding-bottom: 10px;
                }}
                h2 {{
                    color: #34495e;
                    margin-top: 30px;
                }}
                h3 {{
                    color: #7f8c8d;
                }}
                .stats {{
                    background-color: #ecf0f1;
                    padding: 20px;
                    border-radius: 5px;
                    margin: 20px 0;
                }}
                .file-link {{
                    display: inline-block;
                    margin: 5px 10px 5px 0;
                    padding: 8px 15px;
                    background-color: #3498db;
                    color: white;
                    text-decoration: none;
                    border-radius: 5px;
                }}
                .file-link:hover {{
                    background-color: #2980b9;
                }}
                ul {{
                    padding-left: 20px;
                }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>AutoClustal Analysis Report</h1>
                
                <div class="stats">
                    <h2>Summary Statistics</h2>
                    <ul>
                        <li>Total sequences analyzed: {total_sequences}</li>
                        <li>Number of clusters identified: {total_clusters}</li>
                        <li>Average cluster size: {avg_cluster_size:.1f} sequences</li>
                        <li>Largest cluster: {largest_cluster} sequences</li>
                        <li>Average sequence length: {avg_seq_length:.0f} bp/aa</li>
                    </ul>
                </div>
                
                {organism_stats}
                
                <h2>Output Files</h2>
                <p>The following files have been generated:</p>
                
                <h3>Data Tables</h3>
                <a href="summary_table.csv" class="file-link">Summary Table (CSV)</a>
                <a href="summary_table.html" class="file-link">Summary Table (HTML)</a>
                <a href="cluster_statistics.csv" class="file-link">Cluster Statistics</a>
                
                <h3>Sequence Analysis</h3>
                <a href="alignment.fasta" class="file-link">Multiple Sequence Alignment</a>
                <a href="tree.newick" class="file-link">Phylogenetic Tree</a>
                
                <h3>Visualizations</h3>
                <a href="phylogenetic_tree.png" class="file-link">Tree Plot</a>
                <a href="cluster_analysis.png" class="file-link">Cluster Analysis</a>
                <a href="pca_analysis.png" class="file-link">PCA Analysis</a>
                <a href="pca_analysis.html" class="file-link">Interactive PCA</a>
                <a href="annotation_analysis.png" class="file-link">Annotation Analysis</a>
                
                <h2>Analysis Pipeline</h2>
                <ol>
                    <li>Sequence loading and validation</li>
                    <li>Multiple sequence alignment</li>
                    <li>Phylogenetic tree construction</li>
                    <li>Sequence clustering</li>
                    <li>Representative sequence selection</li>
                    <li>Database searches (BLAST/BLAT)</li>
                    <li>PCA and statistical analysis</li>
                    <li>Report generation</li>
                </ol>
                
                <h2>Methodology</h2>
                <p>This analysis was performed using the AutoClustal pipeline, which combines multiple bioinformatics tools and methods for comprehensive sequence analysis. The pipeline includes sequence alignment, phylogenetic reconstruction, clustering analysis, and functional annotation through database searches.</p>
                
                <p><em>Report generated on {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</em></p>
            </div>
        </body>
        </html>
        """
        
        with open(report_path, 'w') as f:
            f.write(html_content)
        
        self.logger.info(f"Final report saved: {report_path}")