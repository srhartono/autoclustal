"""
Simple Analysis Module
=====================

Simplified analysis and reporting module for AutoClustal.
"""

import os
import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Optional, Any

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False

class AnalysisReporter:
    """
    Simple analysis and reporting for sequence analysis results.
    """
    
    def __init__(self, output_dir: Path):
        """Initialize analysis reporter."""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
    
    def generate_summary_table(self, sequences: Dict[str, Any], 
                             clusters: Dict[str, List[str]], 
                             annotations: Dict[str, Any] = None):
        """Generate comprehensive summary table."""
        self.logger.info("Generating summary table")
        
        summary_data = []
        
        for cluster_id, seq_ids in clusters.items():
            for seq_id in seq_ids:
                if seq_id in sequences:
                    seq_record = sequences[seq_id]
                    
                    row = {
                        'sequence_id': seq_id,
                        'cluster_id': cluster_id,
                        'sequence_length': len(str(seq_record.seq).replace('-', '')),
                        'description': getattr(seq_record, 'description', ''),
                    }
                    
                    if annotations and seq_id in annotations:
                        annotation = annotations[seq_id]
                        if 'hits' in annotation and annotation['hits']:
                            best_hit = annotation['hits'][0]
                            row.update({
                                'best_match_title': best_hit.get('title', ''),
                                'best_match_evalue': best_hit.get('evalue', ''),
                                'organism': self._extract_organism_from_hit(best_hit)
                            })
                        else:
                            row.update({
                                'best_match_title': 'No significant hits',
                                'organism': 'Unknown'
                            })
                    else:
                        row.update({
                            'best_match_title': 'Not searched',
                            'organism': 'Unknown'
                        })
                    
                    summary_data.append(row)
        
        # Save summary table
        if PANDAS_AVAILABLE:
            df = pd.DataFrame(summary_data)
            csv_path = self.output_dir / "summary_table.csv"
            df.to_csv(csv_path, index=False)
            self.logger.info(f"Summary table saved: {csv_path}")
        else:
            self._save_basic_table(summary_data, 'summary_table.txt')
    
    def _extract_organism_from_hit(self, hit: Dict[str, Any]) -> str:
        """Extract organism name from BLAST hit."""
        title = hit.get('title', '')
        
        import re
        bracket_match = re.search(r'\[([^\]]+)\]', title)
        if bracket_match:
            return bracket_match.group(1)
        
        return 'Unknown'
    
    def _save_basic_table(self, data: List[Dict], filename: str):
        """Save table in basic text format."""
        filepath = self.output_dir / filename
        
        if not data:
            with open(filepath, 'w') as f:
                f.write("No data to report\n")
            return
        
        headers = sorted(set().union(*(d.keys() for d in data)))
        
        with open(filepath, 'w') as f:
            f.write('\t'.join(headers) + '\n')
            for row in data:
                values = [str(row.get(key, '')) for key in headers]
                f.write('\t'.join(values) + '\n')
        
        self.logger.info(f"Summary table saved: {filepath}")
    
    def perform_pca_analysis(self, alignment, clusters: Dict[str, List[str]], 
                           annotations: Dict[str, Any] = None):
        """Perform basic PCA analysis (placeholder)."""
        self.logger.info("PCA analysis requested - would require scikit-learn")
        
        # Create a simple placeholder plot
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, 'PCA Analysis\n(requires scikit-learn)', 
                ha='center', va='center', fontsize=16)
        plt.axis('off')
        
        pca_path = self.output_dir / 'pca_analysis.png'
        plt.savefig(pca_path)
        plt.close()
        
        self.logger.info(f"PCA placeholder saved: {pca_path}")
    
    def create_visualizations(self, tree=None, clusters: Dict[str, List[str]] = None, 
                            annotations: Dict[str, Any] = None):
        """Create basic visualizations."""
        self.logger.info("Creating visualizations")
        
        if clusters:
            self._create_cluster_plot(clusters)
    
    def _create_cluster_plot(self, clusters: Dict[str, List[str]]):
        """Create cluster size visualization."""
        cluster_sizes = [len(seq_ids) for seq_ids in clusters.values()]
        cluster_names = list(clusters.keys())
        
        plt.figure(figsize=(10, 6))
        plt.bar(cluster_names, cluster_sizes)
        plt.xlabel('Cluster ID')
        plt.ylabel('Number of Sequences')
        plt.title('Cluster Size Distribution')
        plt.xticks(rotation=45)
        plt.tight_layout()
        
        plot_path = self.output_dir / 'cluster_analysis.png'
        plt.savefig(plot_path)
        plt.close()
        
        self.logger.info(f"Cluster plot saved: {plot_path}")