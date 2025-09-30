#!/usr/bin/env python3
"""
Example Usage Script for AutoClustal
===================================

This script demonstrates how to use the AutoClustal pipeline
with example sequences.
"""

import os
import sys
from pathlib import Path

# Add the AutoClustal directory to Python path
autoclustal_dir = Path(__file__).parent
sys.path.insert(0, str(autoclustal_dir))

def create_example_sequences():
    """Create example FASTA file for testing."""
    example_sequences = """
>seq1_human
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC

>seq2_mouse
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC

>seq3_rat
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCAC

>seq4_dog
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC

>seq5_fish
ATGAAAGACATCAAACAACTGGACATGCTGGAGGCCATCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC

>seq6_chicken
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC

>seq7_frog
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGTGC

>seq8_plant
ATGAAAGACGTCAAACAACTGGACATGCTGGAGGCCGTCAAGGAGCTGCTGAAGGAGGCCGAG
GAGGAGCTGGAGCGCCTGGAGAAGGAGCTGGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCT
GGACGAGCTGGAGCGCGAGCTGGAGAAGGAGCTGGACGAGCTGGAGCGC
""".strip()
    
    with open("example_sequences.fasta", "w") as f:
        f.write(example_sequences)
    
    print("Created example_sequences.fasta")

def run_example():
    """Run AutoClustal with example data."""
    print("AutoClustal Example Run")
    print("=" * 50)
    
    # Create example sequences
    create_example_sequences()
    
    # Import AutoClustal modules
    try:
        from modules.sequence_handler import SequenceHandler
        from modules.alignment import AlignmentEngine
        from modules.phylogeny import PhylogeneticAnalyzer
        from modules.clustering import SequenceClusterer
        
        print("✓ All modules imported successfully")
        
        # Initialize components
        seq_handler = SequenceHandler()
        aligner = AlignmentEngine(tool='muscle')
        phylo_analyzer = PhylogeneticAnalyzer(method='nj')
        clusterer = SequenceClusterer(threshold=0.8)
        
        print("✓ Components initialized")
        
        # Load sequences
        print("\n1. Loading sequences...")
        sequences = seq_handler.load_sequences(["example_sequences.fasta"], seq_type='auto')
        print(f"   Loaded {len(sequences)} sequences")
        
        # Get sequence stats
        stats = seq_handler.get_sequence_stats(sequences)
        print(f"   Average length: {stats.get('mean_length', 0):.0f} bp")
        
        # Perform alignment
        print("\n2. Performing alignment...")
        alignment = aligner.align_sequences(sequences, "example_alignment.fasta")
        print("   ✓ Alignment completed")
        
        # Build phylogenetic tree
        print("\n3. Building phylogenetic tree...")
        tree = phylo_analyzer.build_tree(alignment, "example_tree.newick")
        print("   ✓ Tree constructed")
        
        # Cluster sequences
        print("\n4. Clustering sequences...")
        clusters = clusterer.cluster_sequences(alignment, tree)
        print(f"   Found {len(clusters)} clusters")
        
        for cluster_id, seq_ids in clusters.items():
            print(f"   {cluster_id}: {len(seq_ids)} sequences")
        
        # Select representatives
        representatives = clusterer.select_representatives(clusters, sequences)
        print(f"   Selected {len(representatives)} representatives")
        
        print("\n✓ Example pipeline completed successfully!")
        print("\nOutput files:")
        print("  - example_alignment.fasta")
        print("  - example_tree.newick")
        
    except ImportError as e:
        print(f"✗ Import error: {e}")
        print("Make sure all required packages are installed:")
        print("pip install -r requirements.txt")
    except Exception as e:
        print(f"✗ Error during execution: {e}")

if __name__ == "__main__":
    run_example()