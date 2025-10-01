#!/usr/bin/env python3
"""
AutoClustal: Comprehensive Sequence Analysis Pipeline
=====================================================

A bioinformatics tool for                 # Propagate annotations to all sequences in clusters
                logger.info("Propagating annotations to all sequences in clusters...")
                logger.debug(f"Representative annotations: {rep_annotations}")
                for cluster_id, seq_ids in clusters.items():uence alignment, phylogenetic analysis,
and functional annotation of FASTA/FASTQ sequences.

Features:
- Automatic sequence type detection (DNA/RNA/protein)
- Multiple sequence alignment
- Phylogenetic tree construction
- Sequence clustering and representative selection
- BLAST/BLAT database searches
- PCA analysis and visualization
- Comprehensive reporting

Author: Auto-generated bioinformatics pipeline
Date: September 30, 2025
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Add local modules to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from modules.sequence_handler import SequenceHandler
from modules.alignment import AlignmentEngine
from modules.phylogeny import PhylogeneticAnalyzer
from modules.blast_search import BlastSearcher
from modules.clustering import SequenceClusterer
from modules.simple_analysis import AnalysisReporter

def setup_logging(log_level='INFO'):
    """Setup logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('autoclustal.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )

def main():
    """Main pipeline execution."""
    parser = argparse.ArgumentParser(
        description='AutoClustal: Comprehensive Sequence Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python autoclustal.py -i sequences.fasta -o results/
  python autoclustal.py -i sequences.fastq -t protein -a muscle -p ml
  python autoclustal.py -i *.fasta -o analysis/ --blast --pca
        """
    )
    
    # Input/Output arguments
    parser.add_argument('-i', '--input', required=True, nargs='+',
                       help='Input FASTA/FASTQ files (supports wildcards)')
    parser.add_argument('-o', '--output', default='autoclustal_results',
                       help='Output directory (default: autoclustal_results)')
    
    # Sequence type arguments
    parser.add_argument('-t', '--type', choices=['auto', 'dna', 'rna', 'protein'],
                       default='auto', help='Sequence type (default: auto-detect)')
    
    # Alignment arguments
    parser.add_argument('-a', '--aligner', choices=['muscle', 'clustalw', 'mafft'],
                       default='muscle', help='Alignment tool (default: muscle)')
    
    # Phylogeny arguments
    parser.add_argument('-p', '--phylogeny', choices=['nj', 'upgma', 'ml'],
                       default='nj', help='Phylogenetic method (default: neighbor-joining)')
    
    # Analysis arguments
    parser.add_argument('--blast', action='store_true',
                       help='Perform BLAST searches for functional annotation')
    parser.add_argument('--blat', action='store_true',
                       help='Perform BLAT searches (requires local BLAT setup)')
    parser.add_argument('--blast-all', action='store_true',
                       help='BLAST search all sequences (not just representatives)')
    parser.add_argument('--pca', action='store_true',
                       help='Perform PCA analysis')
    parser.add_argument('--cluster-threshold', type=float, default=0.8,
                       help='Clustering threshold for representative selection (default: 0.8)')
    parser.add_argument('--force', action='store_true',
                       help='Force new BLAST searches even if cached results exist')
    parser.add_argument('--batch-blast', action='store_true', default=True,
                       help='Use batch BLAST search for multiple representatives (default: True)')
    parser.add_argument('--no-batch-blast', dest='batch_blast', action='store_false',
                       help='Disable batch BLAST search and use individual searches')
    
    # Other arguments
    parser.add_argument('--threads', type=int, default=4,
                       help='Number of threads to use (default: 4)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level (default: INFO)')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    logger = logging.getLogger(__name__)
    
    try:
        # Create output directory with results structure
        if args.output and not args.output.startswith('results/'):
            final_output_dir = Path('results') / args.output
        else:
            final_output_dir = Path(args.output)
        
        final_output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("Starting AutoClustal pipeline...")
        logger.info(f"Input files: {args.input}")
        logger.info(f"Output directory: {final_output_dir}")
        
        # Initialize pipeline components
        seq_handler = SequenceHandler()
        aligner = AlignmentEngine(tool=args.aligner, threads=args.threads)
        phylo_analyzer = PhylogeneticAnalyzer(method=args.phylogeny)
        clusterer = SequenceClusterer(threshold=args.cluster_threshold)
        reporter = AnalysisReporter(output_dir=final_output_dir)
        
        # Step 1: Load and process sequences
        logger.info("Step 1: Loading and processing sequences...")
        sequences = seq_handler.load_sequences(args.input, seq_type=args.type)
        logger.info(f"Loaded {len(sequences)} sequences")
        
        # Step 2: Sequence alignment
        logger.info("Step 2: Performing multiple sequence alignment...")
        alignment = aligner.align_sequences(sequences, final_output_dir / "alignment.fasta")
        
        # Step 3: Phylogenetic analysis
        logger.info("Step 3: Constructing phylogenetic tree...")
        tree = phylo_analyzer.build_tree(alignment, final_output_dir / "tree.newick")
        
        # Step 4: Sequence clustering
        logger.info("Step 4: Clustering sequences and selecting representatives...")
        clusters = clusterer.cluster_sequences(alignment, tree)
        representatives = clusterer.select_representatives(clusters, sequences)
        
        # Step 5: Functional annotation (optional)
        annotations = {}
        if args.blast or args.blat or args.blast_all:
            logger.info("Step 5: Performing database searches...")
            searcher = BlastSearcher(output_dir=str(final_output_dir), force=args.force)
            
            if args.blast_all:
                # Search ALL sequences (slower but more accurate)
                logger.info(f"BLAST searching all {len(sequences)} sequences...")
                for seq_id, seq_record in sequences.items():
                    if args.blast or args.blast_all:
                        annotations[seq_id] = searcher.blast_search(seq_record)
                    elif args.blat:
                        annotations[seq_id] = searcher.blat_search(seq_record)
            else:
                # Search representative sequences only (faster)
                logger.info(f"BLAST searching {len(representatives)} representative sequences...")
                rep_annotations = {}
                
                if args.blast and len(representatives) > 1 and args.batch_blast:
                    # Use batch BLAST for multiple representatives (more efficient)
                    logger.info("Using batch BLAST search for multiple representatives...")
                    try:
                        # Convert representatives to seq_id -> SeqRecord mapping for batch search
                        rep_seq_dict = {rep_seq.id: rep_seq for rep_seq in representatives.values()}
                        batch_results = searcher.batch_blast_search(rep_seq_dict)
                        
                        # Map batch results to cluster IDs
                        for cluster_id, rep_seq in representatives.items():
                            rep_id = rep_seq.id
                            if rep_id in batch_results:
                                rep_annotations[cluster_id] = batch_results[rep_id]
                                logger.debug(f"Batch BLAST result for {cluster_id}: {batch_results[rep_id]}")
                            else:
                                logger.warning(f"No batch result found for {cluster_id} ({rep_id})")
                                # Fallback to individual search
                                rep_annotations[cluster_id] = searcher.blast_search(rep_seq)
                    except Exception as e:
                        logger.warning(f"Batch BLAST failed: {e}. Falling back to individual searches...")
                        # Fallback to individual searches
                        for cluster_id, rep_seq in representatives.items():
                            rep_annotations[cluster_id] = searcher.blast_search(rep_seq)
                else:
                    # Single representative, batch disabled, or non-BLAST search - use individual searches
                    for cluster_id, rep_seq in representatives.items():
                        if args.blast:
                            result = searcher.blast_search(rep_seq)
                            logger.debug(f"BLAST result for {cluster_id}: {result}")
                            rep_annotations[cluster_id] = result
                        elif args.blat:
                            rep_annotations[cluster_id] = searcher.blat_search(rep_seq)
                
                # Propagate annotations to all sequences in each cluster
                logger.info("Propagating annotations to all sequences in clusters...")
                for cluster_id, seq_ids in clusters.items():
                    # Check if we have annotation for this cluster
                    if cluster_id in rep_annotations:
                        rep_annotation = rep_annotations[cluster_id]
                        logger.debug(f"Propagating annotation from {cluster_id} to {len(seq_ids)} sequences: {rep_annotation}")
                        for seq_id in seq_ids:
                            # Copy the representative's annotation to each sequence in cluster
                            annotations[seq_id] = rep_annotation.copy()
                            # Mark that this was propagated from representative
                            if 'hits' in annotations[seq_id] and annotations[seq_id]['hits']:
                                annotations[seq_id]['search_note'] = f"Annotation from cluster representative: {cluster_id}"
                    else:
                        # No annotation available for this cluster
                        logger.debug(f"No annotation available for {cluster_id}")
                        for seq_id in seq_ids:
                            annotations[seq_id] = {'hits': [], 'search_note': 'No representative annotation available'}
        
        # Step 6: Analysis and reporting
        logger.info("Step 6: Generating analysis reports...")
        reporter.generate_summary_table(sequences, clusters, annotations)
        
        if args.pca:
            reporter.perform_pca_analysis(alignment, clusters, annotations)
        
        reporter.create_visualizations(tree, clusters, annotations)
        
        logger.info("Pipeline completed successfully!")
        logger.info(f"Results saved to: {final_output_dir}")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()