#!/usr/bin/env python3
"""
AutoClustal Tree and Organism Visualizer
========================================

Creates comprehensive visualizations including:
1. Phylogenetic tree visualization
2. Organism distribution bar plots
3. Model organism analysis
4. Summary statistics
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import re
from typing import Dict, List, Tuple, Optional

# Set up matplotlib style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

class AutoClustalVisualizer:
    """
    Comprehensive visualizer for AutoClustal results.
    """
    
    def __init__(self, results_dir: str):
        """Initialize visualizer with results directory."""
        self.results_dir = Path(results_dir)
        self.output_dir = self.results_dir / "visualizations"
        self.output_dir.mkdir(exist_ok=True)
        
        # Model organism categories
        self.organism_categories = {
            'Human': ['homo sapiens', 'human', 'h. sapiens', 'homo_sapiens'],
            'Mouse': ['mus musculus', 'mouse', 'm. musculus', 'mus_musculus'],
            'Rat': ['rattus norvegicus', 'rat', 'r. norvegicus', 'rattus_norvegicus'],
            'C. elegans': ['caenorhabditis elegans', 'c. elegans', 'c_elegans', 'nematode'],
            'Drosophila': ['drosophila melanogaster', 'd. melanogaster', 'fruit fly', 'drosophila'],
            'Zebrafish': ['danio rerio', 'zebrafish', 'd. rerio', 'danio_rerio'],
            'Arabidopsis': ['arabidopsis thaliana', 'a. thaliana', 'arabidopsis'],
            'Yeast': ['saccharomyces cerevisiae', 's. cerevisiae', 'yeast', 'baker yeast'],
            'E. coli': ['escherichia coli', 'e. coli', 'e_coli'],
            'Bacteria': ['bacteria', 'bacterial', 'bacillus', 'streptococcus', 'staphylococcus'],
            'Virus': ['virus', 'viral', 'phage', 'coronavirus', 'influenza'],
            'Fungi': ['fungi', 'fungal', 'candida', 'aspergillus', 'neurospora'],
            'Plant': ['plant', 'plantae', 'viridiplantae', 'rice', 'wheat', 'maize'],
            'Other': []  # Will be filled with remaining organisms
        }
    
    def load_data(self):
        """Load all necessary data files."""
        # Load summary table
        summary_file = self.results_dir / "summary_table.csv"
        if summary_file.exists():
            self.summary_df = pd.read_csv(summary_file)
            print(f"✓ Loaded summary table: {len(self.summary_df)} sequences")
        else:
            raise FileNotFoundError(f"Summary table not found: {summary_file}")
        
        # Load tree file
        tree_file = self.results_dir / "tree.newick"
        if tree_file.exists():
            with open(tree_file, 'r') as f:
                self.tree_newick = f.read().strip()
            print(f"✓ Loaded phylogenetic tree")
        else:
            print("⚠ Tree file not found, skipping tree visualization")
            self.tree_newick = None
    
    def categorize_organisms(self) -> Dict[str, List[str]]:
        """Categorize organisms into model organism groups."""
        organism_counts = {}
        
        # Initialize categories
        for category in self.organism_categories:
            organism_counts[category] = []
        
        # Process each sequence
        for idx, row in self.summary_df.iterrows():
            organism = str(row.get('organism', 'Unknown')).lower().strip()
            
            # If no organism data or unknown, try to infer from sequence ID (demo mode)
            if organism in ['unknown', 'not searched', 'no significant hits', '', 'nan']:
                seq_id = str(row.get('sequence_id', '')).lower()
                organism = self._infer_organism_from_id(seq_id)
            
            if organism in ['unknown', 'not searched', 'no significant hits', '']:
                continue
            
            # Find category for this organism
            categorized = False
            for category, keywords in self.organism_categories.items():
                if category == 'Other':
                    continue
                
                for keyword in keywords:
                    if keyword in organism:
                        organism_counts[category].append(organism)
                        categorized = True
                        break
                
                if categorized:
                    break
            
            # Add to 'Other' if not categorized
            if not categorized:
                organism_counts['Other'].append(organism)
        
        return organism_counts
    
    def _infer_organism_from_id(self, seq_id: str) -> str:
        """Infer organism from sequence ID for demo purposes."""
        seq_id = seq_id.lower()
        
        # Demo organism inference based on sequence names
        if 'human' in seq_id:
            return 'homo sapiens'
        elif 'mouse' in seq_id:
            return 'mus musculus'
        elif 'rat' in seq_id:
            return 'rattus norvegicus'
        elif 'dog' in seq_id:
            return 'canis familiaris'
        elif 'fish' in seq_id:
            return 'danio rerio'
        elif 'chicken' in seq_id:
            return 'gallus gallus'
        elif 'frog' in seq_id:
            return 'xenopus laevis'
        elif 'plant' in seq_id:
            return 'arabidopsis thaliana'
        elif 'yeast' in seq_id:
            return 'saccharomyces cerevisiae'
        elif 'ecoli' in seq_id or 'coli' in seq_id:
            return 'escherichia coli'
        else:
            return 'unknown'
    
    def create_organism_barplot(self):
        """Create comprehensive organism distribution bar plots."""
        print("Creating organism distribution plots...")
        
        # Categorize organisms
        organism_counts = self.categorize_organisms()
        
        # Count sequences per category
        category_counts = {}
        for category, organisms in organism_counts.items():
            category_counts[category] = len(organisms)
        
        # Remove empty categories
        category_counts = {k: v for k, v in category_counts.items() if v > 0}
        
        if not category_counts:
            print("⚠ No organism data found for plotting")
            return
        
        # Create figure with subplots
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('AutoClustal Organism Analysis', fontsize=16, fontweight='bold')
        
        # 1. Model Organisms Bar Plot
        categories = list(category_counts.keys())
        counts = list(category_counts.values())
        
        # Sort by count (descending)
        sorted_data = sorted(zip(categories, counts), key=lambda x: x[1], reverse=True)
        categories_sorted, counts_sorted = zip(*sorted_data) if sorted_data else ([], [])
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(categories_sorted)))
        bars1 = ax1.bar(categories_sorted, counts_sorted, color=colors, alpha=0.8, edgecolor='black')
        ax1.set_title('Sequences by Model Organism Category', fontweight='bold')
        ax1.set_xlabel('Organism Category')
        ax1.set_ylabel('Number of Sequences')
        ax1.tick_params(axis='x', rotation=45)
        
        # Add value labels on bars
        for bar in bars1:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{int(height)}', ha='center', va='bottom', fontweight='bold')
        
        # 2. Pie Chart of Organism Distribution
        if len(categories_sorted) > 1:
            wedges, texts, autotexts = ax2.pie(counts_sorted, labels=categories_sorted, autopct='%1.1f%%',
                                              colors=colors, startangle=90)
            ax2.set_title('Organism Distribution (Percentage)', fontweight='bold')
            
            # Make percentage text bold
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
        else:
            ax2.text(0.5, 0.5, 'Single Category\nFound', ha='center', va='center',
                    transform=ax2.transAxes, fontsize=12)
            ax2.set_title('Organism Distribution (Percentage)', fontweight='bold')
        
        # 3. Top Individual Organisms
        all_organisms = []
        for organisms in organism_counts.values():
            all_organisms.extend(organisms)
        
        if all_organisms:
            organism_individual_counts = {}
            for org in all_organisms:
                organism_individual_counts[org] = organism_individual_counts.get(org, 0) + 1
            
            # Get top 15 individual organisms
            top_organisms = sorted(organism_individual_counts.items(), 
                                 key=lambda x: x[1], reverse=True)[:15]
            
            if top_organisms:
                org_names, org_counts = zip(*top_organisms)
                
                bars3 = ax3.barh(range(len(org_names)), org_counts, 
                               color=plt.cm.viridis(np.linspace(0, 1, len(org_names))))
                ax3.set_yticks(range(len(org_names)))
                ax3.set_yticklabels([name.title() for name in org_names])
                ax3.set_xlabel('Number of Sequences')
                ax3.set_title('Top Individual Organisms', fontweight='bold')
                ax3.invert_yaxis()
                
                # Add value labels
                for i, (bar, count) in enumerate(zip(bars3, org_counts)):
                    ax3.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                            f'{count}', ha='left', va='center', fontweight='bold')
        
        # 4. Summary Statistics
        total_sequences = len(self.summary_df)
        annotated_sequences = len([org for orgs in organism_counts.values() for org in orgs])
        unknown_sequences = total_sequences - annotated_sequences
        
        summary_data = [
            ['Total Sequences', total_sequences],
            ['Annotated Sequences', annotated_sequences],
            ['Unknown/Unannotated', unknown_sequences],
            ['Unique Organisms', len(set(all_organisms))],
            ['Organism Categories', len([cat for cat, count in category_counts.items() if count > 0])]
        ]
        
        ax4.axis('tight')
        ax4.axis('off')
        table = ax4.table(cellText=summary_data,
                         colLabels=['Metric', 'Count'],
                         cellLoc='center',
                         loc='center',
                         colWidths=[0.6, 0.3])
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 1.5)
        
        # Style the table
        for i in range(len(summary_data) + 1):
            for j in range(2):
                cell = table[(i, j)]
                if i == 0:  # Header
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
                    if j == 1:  # Count column
                        cell.set_text_props(weight='bold')
        
        ax4.set_title('Analysis Summary', fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        output_path = self.output_dir / 'organism_analysis.png'
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"✓ Organism analysis saved: {output_path}")
        
        # Also save as PDF for publication
        pdf_path = self.output_dir / 'organism_analysis.pdf'
        plt.savefig(pdf_path, dpi=300, bbox_inches='tight')
        print(f"✓ PDF version saved: {pdf_path}")
        
        #plt.show()
        plt.close()
        
        return category_counts
    
    def create_tree_visualization(self):
        """Create phylogenetic tree visualization."""
        if not self.tree_newick:
            print("⚠ No tree data available")
            return
        
        print("Creating phylogenetic tree visualization...")
        
        try:
            # Try to use BioPython for tree visualization
            from Bio import Phylo
            from io import StringIO
            
            # Parse tree
            tree_handle = StringIO(self.tree_newick)
            tree = Phylo.read(tree_handle, 'newick')
            
            # Create tree plot
            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            
            Phylo.draw(tree, axes=ax, do_show=False)
            ax.set_title('Phylogenetic Tree', fontsize=14, fontweight='bold')
            
            plt.tight_layout()
            
            # Save tree plot
            tree_path = self.output_dir / 'phylogenetic_tree.png'
            plt.savefig(tree_path, dpi=300, bbox_inches='tight')
            print(f"✓ Tree visualization saved: {tree_path}")
            
            #plt.show()
            plt.close()
            
        except ImportError:
            print("⚠ BioPython not available for tree visualization")
            self._create_basic_tree_plot()
        except Exception as e:
            print(f"⚠ Tree visualization failed: {e}")
            self._create_basic_tree_plot()
    
    def _create_basic_tree_plot(self):
        """Create a basic tree representation when BioPython is not available."""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        ax.text(0.5, 0.6, 'Phylogenetic Tree', ha='center', va='center', 
                fontsize=16, fontweight='bold')
        ax.text(0.5, 0.5, f'Newick Format:', ha='center', va='center', 
                fontsize=12, fontweight='bold')
        
        # Display truncated Newick string
        newick_display = self.tree_newick[:200] + "..." if len(self.tree_newick) > 200 else self.tree_newick
        ax.text(0.5, 0.4, newick_display, ha='center', va='center', 
                fontsize=8, family='monospace', wrap=True)
        
        ax.text(0.5, 0.2, '(Install BioPython for enhanced tree visualization)', 
                ha='center', va='center', fontsize=10, style='italic')
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        
        plt.tight_layout()
        
        # Save basic tree plot
        tree_path = self.output_dir / 'phylogenetic_tree_basic.png'
        plt.savefig(tree_path, dpi=300, bbox_inches='tight')
        print(f"✓ Basic tree visualization saved: {tree_path}")
        
        #plt.show()
        plt.close()
    
    def create_cluster_analysis(self):
        """Create cluster analysis visualization."""
        print("Creating cluster analysis...")
        
        if 'cluster_id' not in self.summary_df.columns:
            print("⚠ No cluster information found")
            return
        
        # Cluster analysis
        cluster_counts = self.summary_df['cluster_id'].value_counts()
        
        # Get representative best_match_title for each cluster
        cluster_info = {}
        for cluster_id in cluster_counts.index:
            cluster_data = self.summary_df[self.summary_df['cluster_id'] == cluster_id]
            
            # Get the best match title (first non-empty one or most common one)
            best_matches = cluster_data['best_match_title'].dropna()
            best_matches = best_matches[best_matches != 'No significant hits']
            
            if len(best_matches) > 0:
                # Use the most common best match title for this cluster
                best_match = best_matches.mode().iloc[0] if len(best_matches.mode()) > 0 else best_matches.iloc[0]
                # Truncate long titles
                if len(best_match) > 40:
                    best_match = best_match[:37] + "..."
            else:
                best_match = "No significant hits"
            
            cluster_info[cluster_id] = {
                'count': cluster_counts[cluster_id],
                'best_match': best_match
            }
        
        # Create simple cluster analysis plot
        self._create_simple_cluster_plot(cluster_counts, cluster_info)
        
        # Create detailed cluster analysis plot
        self._create_detailed_cluster_plot(cluster_counts, cluster_info)
    
    def _create_simple_cluster_plot(self, cluster_counts, cluster_info):
        """Create simple cluster analysis plot with best match titles."""
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        # Create labels with cluster ID and best match
        labels = []
        for cluster_id in cluster_counts.index:
            best_match = cluster_info[cluster_id]['best_match']
            labels.append(f"{cluster_id}\n{best_match}")
        
        # Create bar plot
        bars = ax.bar(range(len(cluster_counts)), cluster_counts.values, 
                     color='skyblue', alpha=0.8, edgecolor='black')
        
        ax.set_title('Clusters with Best Match Annotations', fontweight='bold', fontsize=14)
        ax.set_xlabel('Cluster ID + Best Match Title')
        ax.set_ylabel('Number of Sequences')
        ax.set_xticks(range(len(cluster_counts)))
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        
        # Add value labels on bars
        for i, (bar, count) in enumerate(zip(bars, cluster_counts.values)):
            ax.text(bar.get_x() + bar.get_width()/2., count + 0.1,
                   str(count), ha='center', va='bottom', fontweight='bold')
        
        plt.tight_layout()
        
        # Save simple cluster analysis
        cluster_path = self.output_dir / 'cluster_analysis.png'
        plt.savefig(cluster_path, dpi=300, bbox_inches='tight')
        print(f"✓ Simple cluster analysis saved: {cluster_path}")
        
        plt.close()
    
    def _create_detailed_cluster_plot(self, cluster_counts, cluster_info):
        """Create detailed cluster analysis with multiple views."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Detailed Cluster Analysis with Annotations', fontsize=16, fontweight='bold')
        
        # 1. Cluster size distribution with best matches
        labels = []
        for cluster_id in cluster_counts.index:
            best_match = cluster_info[cluster_id]['best_match']
            # Even shorter labels for this plot
            short_match = best_match[:25] + "..." if len(best_match) > 25 else best_match
            labels.append(f"{cluster_id}\n{short_match}")
        
        bars1 = ax1.bar(range(len(cluster_counts)), cluster_counts.values, 
                       color='skyblue', alpha=0.8, edgecolor='black')
        ax1.set_title('Sequences per Cluster (with Best Match)', fontweight='bold')
        ax1.set_xlabel('Cluster ID + Annotation')
        ax1.set_ylabel('Number of Sequences')
        ax1.set_xticks(range(len(cluster_counts)))
        ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=8)
        
        # Add value labels
        for i, (bar, count) in enumerate(zip(bars1, cluster_counts.values)):
            ax1.text(bar.get_x() + bar.get_width()/2., count + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')
        
        # 2. Cluster size histogram
        ax2.hist(cluster_counts.values, bins=max(1, len(cluster_counts)//2), 
                alpha=0.7, color='lightcoral', edgecolor='black')
        ax2.set_title('Cluster Size Distribution', fontweight='bold')
        ax2.set_xlabel('Cluster Size (sequences)')
        ax2.set_ylabel('Number of Clusters')
        
        # 3. Annotation success rate per cluster
        annotated_clusters = 0
        unannotated_clusters = 0
        
        for cluster_id, info in cluster_info.items():
            if info['best_match'] != "No significant hits":
                annotated_clusters += 1
            else:
                unannotated_clusters += 1
        
        if annotated_clusters > 0 or unannotated_clusters > 0:
            sizes = [annotated_clusters, unannotated_clusters]
            labels_pie = ['Annotated Clusters', 'Unannotated Clusters']
            colors = ['lightgreen', 'lightcoral']
            
            wedges, texts, autotexts = ax3.pie(sizes, labels=labels_pie, autopct='%1.1f%%',
                                              colors=colors, startangle=90)
            ax3.set_title('Cluster Annotation Success Rate', fontweight='bold')
            
            # Make percentage text bold
            for autotext in autotexts:
                autotext.set_color('white')
                autotext.set_fontweight('bold')
        
        # 4. Summary table
        total_clusters = len(cluster_counts)
        total_sequences = sum(cluster_counts.values)
        avg_cluster_size = total_sequences / total_clusters if total_clusters > 0 else 0
        largest_cluster = max(cluster_counts.values) if len(cluster_counts) > 0 else 0
        smallest_cluster = min(cluster_counts.values) if len(cluster_counts) > 0 else 0
        
        summary_data = [
            ['Total Clusters', total_clusters],
            ['Total Sequences', total_sequences],
            ['Avg Cluster Size', f'{avg_cluster_size:.1f}'],
            ['Largest Cluster', largest_cluster],
            ['Smallest Cluster', smallest_cluster],
            ['Annotated Clusters', annotated_clusters],
            ['Unannotated Clusters', unannotated_clusters]
        ]
        
        ax4.axis('tight')
        ax4.axis('off')
        table = ax4.table(cellText=summary_data,
                         colLabels=['Metric', 'Value'],
                         cellLoc='center',
                         loc='center',
                         colWidths=[0.6, 0.4])
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.5)
        
        # Style the table
        for i in range(len(summary_data) + 1):
            for j in range(2):
                cell = table[(i, j)]
                if i == 0:  # Header
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
                    if j == 1:  # Value column
                        cell.set_text_props(weight='bold')
        
        ax4.set_title('Cluster Summary Statistics', fontweight='bold')
        
        plt.tight_layout()
        
        # Save detailed cluster analysis
        cluster_path = self.output_dir / 'cluster_analysis_detailed.png'
        plt.savefig(cluster_path, dpi=300, bbox_inches='tight')
        print(f"✓ Detailed cluster analysis saved: {cluster_path}")
        
        plt.close()
    
    def generate_summary_report(self, organism_counts: Dict[str, int]):
        """Generate a comprehensive text summary report."""
        print("Generating summary report...")
        
        total_sequences = len(self.summary_df)
        annotated_sequences = sum(organism_counts.values()) if organism_counts else 0
        
        report_lines = [
            "=" * 60,
            "AUTOCLUSTAL ANALYSIS SUMMARY REPORT",
            "=" * 60,
            f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"Results Directory: {self.results_dir}",
            "",
            "SEQUENCE OVERVIEW",
            "-" * 20,
            f"Total Sequences Analyzed: {total_sequences}",
            f"Successfully Annotated: {annotated_sequences}",
            f"Unknown/Unannotated: {total_sequences - annotated_sequences}",
            f"Annotation Success Rate: {(annotated_sequences/total_sequences)*100:.1f}%" if total_sequences > 0 else "No sequences",
            "",
            "ORGANISM DISTRIBUTION",
            "-" * 20
        ]
        
        # Add organism counts
        if organism_counts:
            for category, count in sorted(organism_counts.items(), key=lambda x: x[1], reverse=True):
                if count > 0:
                    percentage = (count / annotated_sequences) * 100 if annotated_sequences > 0 else 0
                    report_lines.append(f"{category:15}: {count:3d} sequences ({percentage:5.1f}%)")
        else:
            report_lines.append("No organism data available")
        
        # Add cluster information if available
        if 'cluster_id' in self.summary_df.columns:
            cluster_counts = self.summary_df['cluster_id'].value_counts()
            report_lines.extend([
                "",
                "CLUSTERING ANALYSIS",
                "-" * 20,
                f"Number of Clusters: {len(cluster_counts)}",
                f"Largest Cluster: {cluster_counts.max()} sequences",
                f"Average Cluster Size: {cluster_counts.mean():.1f} sequences"
            ])
        
        report_lines.extend([
            "",
            "OUTPUT FILES GENERATED",
            "-" * 20,
            f"• Organism Analysis: {self.output_dir}/organism_analysis.png",
            f"• Phylogenetic Tree: {self.output_dir}/phylogenetic_tree.png",
            f"• Cluster Analysis: {self.output_dir}/cluster_analysis_detailed.png",
            f"• Summary Report: {self.output_dir}/analysis_summary.txt",
            "",
            "=" * 60
        ])
        
        # Save report
        report_path = self.output_dir / 'analysis_summary.txt'
        with open(report_path, 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"✓ Summary report saved: {report_path}")
        
        # Also print to console
        print('\n'.join(report_lines))
    
    def run_visualization(self):
        """Run complete visualization pipeline."""
        print("Starting AutoClustal Visualization Pipeline...")
        print("=" * 50)
        
        # Load data
        self.load_data()
        
        # Create visualizations
        organism_counts = self.create_organism_barplot()
        self.create_tree_visualization()
        self.create_cluster_analysis()
        
        # Generate summary report
        self.generate_summary_report(organism_counts)
        
        print("\n" + "=" * 50)
        print("✓ Visualization pipeline completed successfully!")
        print(f"✓ All outputs saved to: {self.output_dir}")

def main():
    """Main function for command-line usage."""
    parser = argparse.ArgumentParser(
        description='AutoClustal Results Visualizer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python grapher.py results/
  python grapher.py results_blast/
  python grapher.py /path/to/autoclustal/output/
        """
    )
    
    parser.add_argument('results_dir', 
                       help='Path to AutoClustal results directory')
    parser.add_argument('--output', '-o',
                       help='Output directory for visualizations (default: results_dir/visualizations)')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.results_dir):
        print(f"Error: Results directory not found: {args.results_dir}")
        sys.exit(1)
    
    # Create visualizer
    visualizer = AutoClustalVisualizer(args.results_dir)
    
    # Override output directory if specified
    if args.output:
        visualizer.output_dir = Path(args.output)
        visualizer.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run visualization
    try:
        visualizer.run_visualization()
    except Exception as e:
        print(f"Error: Visualization failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()