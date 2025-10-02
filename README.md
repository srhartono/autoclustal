# AutoClustal: Comprehensive Sequence Analysis Pipeline

AutoClustal is a powerful bioinformatics tool that performs comprehensive sequence analysis including multiple sequence alignment, phylogenetic analysis, clustering, and functional annotation.

## Features

- **Sequence Input**: Support for FASTA and FASTQ formats
- **Sequence Sampling**: Random sampling from large datasets (-M and -N parameters)
- **Automatic Detection**: DNA, RNA, or protein sequence type detection
- **Multiple Alignment**: Integration with MUSCLE, ClustalW, and MAFFT
- **Phylogenetic Analysis**: Distance-based (NJ, UPGMA) and ML methods
- **Sequence Clustering**: Multiple clustering algorithms with representative selection
- **Database Searches**: BLAST and BLAT integration for functional annotation
- **Statistical Analysis**: PCA, t-SNE, and comprehensive reporting (WIP)
- **Advanced Visualizations**: 
  - Phylogenetic tree plots
  - Organism distribution bar charts
  - Model organism analysis (Human, Mouse, C. elegans, Yeast, E. coli, etc.)
  - Cluster analysis plots
  - Publication-ready figures (PNG + PDF)
- **Interactive Reports**: Comprehensive HTML reports with organism statistics

## Installation

### Option 1: Python Virtual Environment (Recommended)

1. **Clone or download the AutoClustal project**
2. **Set up the virtual environment**:
   ```bash
   # Linux/Mac
   ./setup_env.sh
   
   # Windows
   setup_env.bat
   ```

3. **Activate the virtual environment**:
   ```bash
   # Linux/Mac
   source .venv/bin/activate
   
   # Windows
   .venv\Scripts\activate.bat
   ```

### Option 2: Conda Environment

```bash
# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate autoclustal
```

### Option 3: Manual Installation

#### Requirements

- Python 3.7+
- Required packages (install with `pip install -r requirements.txt`):
  - biopython>=1.81
  - matplotlib>=3.7.0
  - seaborn>=0.12.0
  - pandas>=2.0.0
  - numpy>=1.24.0
  - scikit-learn>=1.3.0
  - scipy>=1.10.0
  - plotly>=5.15.0

### Nextflow Installation (Optional)

For running the Nextflow pipeline, install Nextflow:

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

# Or on Windows, download from https://www.nextflow.io/
```

### Optional External Tools

For optimal performance, install these external tools:

- **MUSCLE**: Multiple sequence alignment
  - Download from: http://www.drive5.com/muscle/
  - Add to PATH

- **ClustalW**: Alternative alignment tool
  - Download from: http://www.clustal.org/clustal2/
  - Add to PATH

- **MAFFT**: Fast alignment tool
  - Download from: https://mafft.cbrc.jp/alignment/software/
  - Add to PATH

- **BLAST+**: Database searches
  - Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  - Add to PATH

## Quick Start

### Option 1: Python Pipeline

```bash
# Activate virtual environment
source .venv/bin/activate  # Linux/Mac
# or
.venv\Scripts\activate.bat  # Windows

# Run with example sequences
python example.py

# Basic analysis
python autoclustal.py -i sequences.fasta -o results/

# Advanced analysis with BLAST and PCA
python autoclustal.py -i sequences.fasta -o results/ --blast --pca
```

### Option 2: Nextflow Pipeline (WIP)

```bash
# Basic Nextflow run
nextflow run autoclustal.nf --input sequences.fasta --outdir results/

# Advanced analysis
nextflow run autoclustal.nf \
  --input "*.fasta" \
  --outdir results \
  --blast \
  --pca \
  --aligner mafft \
  --phylogeny nj

# Run with specific profile
nextflow run autoclustal.nf --input sequences.fasta -profile docker
```

### Option 3: Visualization Tools

```bash
# Quick visualization of results
python visualize.py results/

# Advanced visualization with custom options
python grapher.py results/ 

# Create demo data for testing
python create_demo_data.py
python visualize.py demo_results/
```

```bash
# Basic Nextflow run
nextflow run autoclustal.nf --input sequences.fasta --outdir results/

# Advanced analysis
nextflow run autoclustal.nf \
  --input "*.fasta" \
  --outdir results \
  --blast \
  --pca \
  --aligner mafft \
  --phylogeny nj

# Run with specific profile
nextflow run autoclustal.nf --input sequences.fasta -profile docker
```

### Command Line Options

#### Python Pipeline
```bash
python autoclustal.py [options]

Required:
  -i, --input           Input FASTA/FASTQ files (supports wildcards)

Optional:
  -o, --output          Output directory (default: autoclustal_results)
  -M, --max-sequences   Maximum number of sequences to take from beginning (default: all)
  -N, --random-sample   Number of sequences to randomly sample from first M (default: all)
  --random-seed         Random seed for reproducible sampling (default: None)
  -t, --type            Sequence type: auto, dna, rna, protein (default: auto)
  -a, --aligner         Alignment tool: muscle, clustalw, mafft (default: muscle)
  -p, --phylogeny       Phylogenetic method: nj, upgma, ml (default: nj)
  --blast               Perform BLAST searches
  --blat                Perform BLAT searches (requires local setup)
  --pca                 Perform PCA analysis
  --cluster-threshold   Clustering threshold (default: 0.8)
  --threads             Number of threads (default: 4)
  --log-level           Logging level: DEBUG, INFO, WARNING, ERROR
```

#### Nextflow Pipeline
```bash
nextflow run autoclustal.nf [options]

Required:
  --input               Input FASTA/FASTQ files

Optional:
  --outdir              Output directory (default: results)
  --seq_type            Sequence type: auto, dna, rna, protein (default: auto)
  --aligner             Alignment tool: muscle, clustalw, mafft (default: muscle)
  --phylogeny           Phylogenetic method: nj, upgma, ml (default: nj)
  --blast               Enable BLAST searches (default: false)
  --blat                Enable BLAT searches (default: false)
  --pca                 Enable PCA analysis (default: false)
  --cluster_threshold   Clustering threshold (default: 0.8)
  --threads             Number of threads (default: 4)
  --max_memory          Maximum memory (default: 16.GB)
  --max_cpus            Maximum CPUs (default: 8)
  --max_time            Maximum time (default: 8.h)
```

### Examples

#### Example 1: Basic DNA sequence analysis
```bash
python autoclustal.py -i dna_sequences.fasta -t dna -a muscle -p nj -o dna_analysis/
```

#### Example 2: Protein analysis with BLAST annotation
```bash
python autoclustal.py -i proteins.fasta -t protein --blast --pca -o protein_analysis/
```

#### Example 3: Multiple files with clustering
```bash
python autoclustal.py -i *.fasta -a mafft --cluster-threshold 0.7 -o multi_analysis/
```

#### Example 4: Large dataset sampling (analyze 100 sequences from first 100,000)
```bash
python autoclustal.py -i huge_dataset.fasta -M 100000 -N 100 --random-seed 42 -o sampled_analysis/
```

#### Example 5: Random sampling from entire dataset
```bash
python autoclustal.py -i large_file.fastq -N 500 --blast --pca -o random_sample_analysis/
```

## Output Files

AutoClustal generates comprehensive output including:

### Data Tables
- `summary_table.csv/html` - Complete analysis summary
- `cluster_statistics.csv` - Cluster-level statistics
- `pca_results.csv` - PCA analysis results (if requested)

### Sequence Analysis
- `alignment.fasta` - Multiple sequence alignment
- `tree.newick` - Phylogenetic tree in Newick format

### Visualizations
- `organism_analysis.png/pdf` - Comprehensive organism distribution plots:
  - Model organism bar chart (Human, Mouse, C. elegans, Yeast, etc.)
  - Organism distribution pie chart
  - Top individual organisms
  - Analysis summary table
- `phylogenetic_tree.png` - Phylogenetic tree visualization
- `cluster_analysis_detailed.png` - Cluster size distributions
- `pca_analysis.png/html` - PCA plots (static and interactive)
- `tsne_analysis.png` - t-SNE visualization
- `annotation_analysis.png` - BLAST/BLAT results summary

### Reports
- `AutoClustal_Report.html` - Comprehensive HTML report
- `analysis_summary.txt` - Text-based summary with organism statistics

## Methodology

### 1. Sequence Processing
- Automatic format detection (FASTA/FASTQ)
- Sequence type detection (DNA/RNA/protein)
- Quality validation and filtering

### 2. Multiple Sequence Alignment
- Support for MUSCLE, ClustalW, and MAFFT
- Automatic parameter optimization
- Alignment quality assessment

### 3. Phylogenetic Analysis
- Distance-based methods (Neighbor-Joining, UPGMA)
- Maximum likelihood (basic implementation)
- Bootstrap support (when available)

### 4. Sequence Clustering
- Distance-based hierarchical clustering
- K-means clustering with optimal k selection
- DBSCAN for density-based clustering
- Representative sequence selection

### 5. Functional Annotation
- BLAST searches against public databases
- Local BLAT searches (when configured)
- Organism identification and classification
- Top hit analysis and statistics

### 6. Statistical Analysis
- Principal Component Analysis (PCA)
- t-SNE for non-linear dimensionality reduction
- Sequence composition analysis
- Clustering validation metrics

## Troubleshooting

### Common Issues

1. **External tools not found**
   - Ensure tools are installed and in PATH
   - Use `which muscle` (Linux/Mac) or `where muscle` (Windows) to verify

2. **BioPython import errors**
   - Install with: `pip install biopython`
   - Update to latest version: `pip install --upgrade biopython`

3. **Memory issues with large datasets**
   - Reduce number of sequences
   - Use sampling for initial analysis
   - Increase system memory

4. **BLAST searches failing**
   - Check internet connection
   - Try local BLAST setup
   - Reduce number of sequences for online searches

### Performance Tips

- Use MAFFT for large alignments (fastest)
- Enable multithreading with `--threads`
- Use representative sequences for initial analysis
- Consider sequence sampling for very large datasets

## API Usage

AutoClustal can also be used as a Python library:

```python
from modules.sequence_handler import SequenceHandler
from modules.alignment import AlignmentEngine
from modules.phylogeny import PhylogeneticAnalyzer
from modules.clustering import SequenceClusterer
from modules.blast_search import BlastSearcher
from modules.analysis import AnalysisReporter

# Initialize components
seq_handler = SequenceHandler()
aligner = AlignmentEngine()
phylo = PhylogeneticAnalyzer()

# Load and process sequences
sequences = seq_handler.load_sequences(['input.fasta'])
alignment = aligner.align_sequences(sequences)
tree = phylo.build_tree(alignment)
```


## License
MIT License

Copyright (C) 2025 Stella R Hartono

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

