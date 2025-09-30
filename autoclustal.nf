#!/usr/bin/env nextflow

/*
========================================================================================
    AutoClustal Nextflow Pipeline
========================================================================================
    Comprehensive sequence analysis pipeline for FASTA/FASTQ sequences
    including alignment, phylogeny, clustering, and functional annotation.
========================================================================================
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

params.input        = null    // Input sequences (FASTA/FASTQ)
params.outdir       = 'results'
params.seq_type     = 'auto'  // auto, dna, rna, protein
params.aligner      = 'muscle' // muscle, clustalw, mafft
params.phylogeny    = 'nj'     // nj, upgma, ml
params.cluster_threshold = 0.8
params.blast        = false
params.blat         = false
params.pca          = false
params.threads      = 4
params.help         = false
params.max_memory   = '8.GB'
params.max_cpus     = 4
params.max_time     = '4.h'

// Validate required parameters
if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input) {
    error "Input sequences not specified! Use --input /path/to/sequences.fasta"
}

/*
========================================================================================
    HELP MESSAGE
========================================================================================
*/

def helpMessage() {
    log.info"""
    ========================================================================================
                              AutoClustal Nextflow Pipeline
    ========================================================================================
    
    Usage:
      nextflow run autoclustal.nf --input sequences.fasta --outdir results
    
    Required Arguments:
      --input               Path to input FASTA/FASTQ file(s)
    
    Optional Arguments:
      --outdir              Output directory (default: results)
      --seq_type            Sequence type: auto, dna, rna, protein (default: auto)
      --aligner             Alignment tool: muscle, clustalw, mafft (default: muscle)
      --phylogeny           Phylogenetic method: nj, upgma, ml (default: nj)
      --cluster_threshold   Clustering threshold (default: 0.8)
      --blast               Enable BLAST searches (default: false)
      --blat                Enable BLAT searches (default: false)
      --pca                 Enable PCA analysis (default: false)
      --threads             Number of threads (default: 4)
      --max_memory          Maximum memory (default: 8.GB)
      --max_cpus            Maximum CPUs (default: 4)
      --max_time            Maximum time (default: 4.h)
    
    Examples:
      nextflow run autoclustal.nf --input sequences.fasta
      nextflow run autoclustal.nf --input "*.fasta" --blast --pca
      nextflow run autoclustal.nf --input proteins.fasta --seq_type protein --aligner mafft
    ========================================================================================
    """.stripIndent()
}

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {
    
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)
    
    // Stage 1: Load and validate sequences
    LOAD_SEQUENCES(input_ch)
    
    // Stage 2: Multiple sequence alignment
    ALIGN_SEQUENCES(LOAD_SEQUENCES.out.sequences)
    
    // Stage 3: Build phylogenetic tree
    BUILD_TREE(ALIGN_SEQUENCES.out.alignment)
    
    // Stage 4: Cluster sequences
    CLUSTER_SEQUENCES(ALIGN_SEQUENCES.out.alignment, BUILD_TREE.out.tree)
    
    // Stage 5: BLAST/BLAT searches (optional)
    if (params.blast || params.blat) {
        SEARCH_DATABASES(CLUSTER_SEQUENCES.out.representatives)
        annotations = SEARCH_DATABASES.out.annotations
    } else {
        annotations = Channel.empty()
    }
    
    // Stage 6: Generate reports and visualizations
    GENERATE_REPORTS(
        LOAD_SEQUENCES.out.sequences,
        ALIGN_SEQUENCES.out.alignment,
        BUILD_TREE.out.tree,
        CLUSTER_SEQUENCES.out.clusters,
        annotations.ifEmpty(file('NO_FILE'))
    )
}

/*
========================================================================================
    PROCESSES
========================================================================================
*/

process LOAD_SEQUENCES {
    tag "$sequences"
    
    input:
    path sequences
    
    output:
    path "sequences.json", emit: sequences
    path "sequence_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.sequence_handler import SequenceHandler
    import json
    import pickle
    
    # Load sequences
    handler = SequenceHandler()
    sequences = handler.load_sequences(['${sequences}'], seq_type='${params.seq_type}')
    
    # Save sequences as pickle for downstream processes
    with open('sequences.pkl', 'wb') as f:
        pickle.dump(sequences, f)
    
    # Save metadata as JSON
    seq_data = {}
    for seq_id, seq_record in sequences.items():
        seq_data[seq_id] = {
            'id': seq_id,
            'length': len(str(seq_record.seq)),
            'description': str(seq_record.description)
        }
    
    with open('sequences.json', 'w') as f:
        json.dump(seq_data, f, indent=2)
    
    # Generate statistics
    stats = handler.get_sequence_stats(sequences)
    with open('sequence_stats.txt', 'w') as f:
        f.write(f"Total sequences: {stats['count']}\\n")
        f.write(f"Min length: {stats['min_length']}\\n")
        f.write(f"Max length: {stats['max_length']}\\n")
        f.write(f"Mean length: {stats['mean_length']:.1f}\\n")
        f.write(f"Total length: {stats['total_length']}\\n")
    """
}

process ALIGN_SEQUENCES {
    tag "alignment"
    publishDir "${params.outdir}/alignment", mode: 'copy'
    
    cpus params.threads
    memory params.max_memory
    time params.max_time
    
    input:
    path sequences_json
    
    output:
    path "alignment.fasta", emit: alignment
    path "alignment_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.sequence_handler import SequenceHandler
    from modules.alignment import AlignmentEngine
    import pickle
    import json
    
    # Load sequences
    with open('sequences.pkl', 'rb') as f:
        sequences = pickle.load(f)
    
    # Perform alignment
    aligner = AlignmentEngine(tool='${params.aligner}', threads=${params.threads})
    alignment = aligner.align_sequences(sequences, 'alignment.fasta')
    
    # Calculate alignment statistics
    stats = aligner.calculate_alignment_stats(alignment)
    with open('alignment_stats.txt', 'w') as f:
        f.write(f"Alignment length: {stats['alignment_length']}\\n")
        f.write(f"Number of sequences: {stats['num_sequences']}\\n")
        f.write(f"Gap percentage: {stats['gap_percentage']:.2f}%\\n")
        f.write(f"Conservation percentage: {stats['conservation_percentage']:.2f}%\\n")
    """
}

process BUILD_TREE {
    tag "phylogeny"
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    
    memory params.max_memory
    time params.max_time
    
    input:
    path alignment
    
    output:
    path "tree.newick", emit: tree
    path "tree_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.phylogeny import PhylogeneticAnalyzer
    from modules.sequence_handler import SequenceHandler
    
    # Load alignment
    handler = SequenceHandler()
    alignment = handler.load_sequences(['${alignment}'])
    
    # Build tree
    phylo = PhylogeneticAnalyzer(method='${params.phylogeny}')
    tree = phylo.build_tree(alignment, 'tree.newick')
    
    # Calculate tree statistics
    stats = phylo.calculate_tree_stats(tree)
    with open('tree_stats.txt', 'w') as f:
        for key, value in stats.items():
            f.write(f"{key}: {value}\\n")
    """
}

process CLUSTER_SEQUENCES {
    tag "clustering"
    publishDir "${params.outdir}/clustering", mode: 'copy'
    
    memory params.max_memory
    time params.max_time
    
    input:
    path alignment
    path tree
    
    output:
    path "clusters.json", emit: clusters
    path "representatives.fasta", emit: representatives
    path "cluster_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.clustering import SequenceClusterer
    from modules.sequence_handler import SequenceHandler
    import json
    import pickle
    
    # Load data
    handler = SequenceHandler()
    alignment = handler.load_sequences(['${alignment}'])
    
    # Load original sequences
    with open('sequences.pkl', 'rb') as f:
        sequences = pickle.load(f)
    
    # Cluster sequences
    clusterer = SequenceClusterer(threshold=${params.cluster_threshold})
    clusters = clusterer.cluster_sequences(alignment)
    
    # Select representatives
    representatives = clusterer.select_representatives(clusters, sequences)
    
    # Save clusters
    with open('clusters.json', 'w') as f:
        json.dump(clusters, f, indent=2)
    
    # Save representatives
    handler.save_sequences(representatives, 'representatives.fasta')
    
    # Calculate cluster statistics
    stats = clusterer.calculate_cluster_stats(clusters, sequences)
    with open('cluster_stats.txt', 'w') as f:
        for cluster_id, cluster_stats in stats.items():
            f.write(f"{cluster_id}:\\n")
            for key, value in cluster_stats.items():
                f.write(f"  {key}: {value}\\n")
    """
}

process SEARCH_DATABASES {
    tag "database_search"
    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    memory params.max_memory
    time params.max_time
    
    input:
    path representatives
    
    output:
    path "annotations.json", emit: annotations
    path "search_summary.txt", emit: summary
    
    when:
    params.blast || params.blat
    
    script:
    def search_type = params.blast ? "blast" : "blat"
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.blast_search import BlastSearcher
    from modules.sequence_handler import SequenceHandler
    import json
    
    # Load representative sequences
    handler = SequenceHandler()
    representatives = handler.load_sequences(['${representatives}'])
    
    # Perform searches
    searcher = BlastSearcher()
    annotations = {}
    
    for seq_id, seq_record in representatives.items():
        if '${search_type}' == 'blast':
            result = searcher.blast_search(seq_record)
        else:
            result = searcher.blat_search(seq_record)
        annotations[seq_id] = result
    
    # Save annotations
    with open('annotations.json', 'w') as f:
        json.dump(annotations, f, indent=2, default=str)
    
    # Generate summary
    summary = searcher.summarize_search_results(list(annotations.values()))
    with open('search_summary.txt', 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\\n")
    """
}

process GENERATE_REPORTS {
    tag "reporting"
    publishDir "${params.outdir}", mode: 'copy'
    
    memory params.max_memory
    time params.max_time
    
    input:
    path sequences_json
    path alignment
    path tree
    path clusters_json
    path annotations
    
    output:
    path "summary_table.*", emit: tables
    path "*.png", emit: plots
    path "*.html", emit: reports
    
    script:
    def pca_flag = params.pca ? "--pca" : ""
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append('${projectDir}')
    from modules.simple_analysis import AnalysisReporter
    from modules.sequence_handler import SequenceHandler
    import json
    import pickle
    from pathlib import Path
    
    # Load data
    with open('sequences.pkl', 'rb') as f:
        sequences = pickle.load(f)
    
    handler = SequenceHandler()
    alignment = handler.load_sequences(['${alignment}'])
    
    with open('${clusters_json}', 'r') as f:
        clusters = json.load(f)
    
    # Load annotations if available
    annotations = None
    if '${annotations}' != 'NO_FILE':
        with open('${annotations}', 'r') as f:
            annotations = json.load(f)
    
    # Generate reports
    reporter = AnalysisReporter(Path('.'))
    reporter.generate_summary_table(sequences, clusters, annotations)
    
    if '${pca_flag}':
        reporter.perform_pca_analysis(alignment, clusters, annotations)
    
    reporter.create_visualizations(tree=None, clusters=clusters, annotations=annotations)
    """
}

/*
========================================================================================
    COMPLETION
========================================================================================
*/

workflow.onComplete {
    log.info """
    ========================================================================================
                              AutoClustal Pipeline Complete!
    ========================================================================================
    Results directory: ${params.outdir}
    
    Output files:
    - Sequences: ${params.outdir}/sequences/
    - Alignment: ${params.outdir}/alignment/
    - Phylogeny: ${params.outdir}/phylogeny/
    - Clustering: ${params.outdir}/clustering/
    - Annotations: ${params.outdir}/annotations/ (if requested)
    - Reports: ${params.outdir}/
    
    Pipeline completed at: ${new Date()}
    ========================================================================================
    """.stripIndent()
}