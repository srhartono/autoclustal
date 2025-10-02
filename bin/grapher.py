#!/usr/bin/env python3
"""
AutoClustal Tree and Organism Visualizer
========================================

Creates comprehensive visualizations including:
1. Phylogenetic tree visualization
2. Organism distribution bar plots
3. Model organism analysis
4. Summary statistics

Created by Stella R. Hartono 2025

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
from collections import Counter

# Set up matplotlib style
plt.style.use('default')
sns.set_palette("husl")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10


# ---- Compact regexes by category (case-insensitive) ----
GROUP_PATTERNS = {
    "human": r"(?:human|homo\s? sapiens|h. sapiens|homo_sapiens)",
    "mouse": r"(?:mus\s?musculus|mouse|m. musculus|mus_musculus)",
    "rat": r"(?:rattus\s?norvegicus|rat|r. norvegicus|rattus_norvegicus)",
    "zebrafish": r"(?:danio\s?rerio|zebrafish|d. rerio|danio_rerio)",
    "celegans": r"(?:caenorhabditis\s?elegans|c\.?\s?elegans|c_elegans|nematode)",
    "drosophila": r"(?:drosophila\s?melanogaster|d\.?\s?melanogaster|fruit fly|drosophila)",
    "arabidopsis": r"(?:arabidopsis\s?thaliana|a\.?\s?thaliana|arabidopsis)",
    "yeast": r"(?:saccharomyces\s?cerevisiae|s\.?\s?cerevisiae|yeast|baker yeast|saccharomyces)",
    "chicken": r"(?:gallus\s?gallus|chicken|g. gallus|gallus_gallus)",
    "frog": r"(?:xenopus\s?laevis|frog|x. laevis|xenopus_laevis)",
    "ecoli": r"(?:escherichia\s?coli|e\.?\s?coli|e_coli|ecoli)",

    "mammals": r"(?:homo[_\s]?sapiens|h\.?\s?sapiens|human|mus\s?musculus|mouse|murine|"
               r"rattus\s?norvegicus|rat|canis|dog|felis|cat|bos|cow|sus|pig|equus|horse|"
               r"ovis|sheep|capra|goat|oryctolagus|rabbit|macaque|macaca|chimp|"
               r"pan\s?troglodytes|gorilla|pongo|orangutan|dolphin|delphin|whale|balaen|"
               r"elephant|loxodonta|monkey)",
    "fish": r"(?:fish|fishes|ichthy|teleost|salmon|trout|tuna|zebrafish|danio rerio)",
    "invertebrate": r"(?:invertebrate|insect|arthropod|annelid|mollusk|cnidarian|"
                    r"nematode|echinoderm|porifera|platyhelminth|rotifer|tardigrade)",
    "bacteria": r"(?:bacteria|bacterial|bacillus|streptococcus|staphylococcus|lactobacillus|"
                      r"clostridium|corynebacterium|enterobacter|serratia|klebsiella|"
                      r"pseudomonas|salmonella|shigella|vibrio|neisseria|haemophilus|"
                      r"legionella|mycobacterium|listeria|bordetella:bacter|philus |zobium |coccus|cillus|ella|zoogloea |myco|spir|archae|therm|methanogen|halobacterium|"
                      r"pseudomonas|escherichia|shigella|vibrio|campylobacter|onia |bacter|proteus|neisseria|clostridium)",
    "virus": r"(?:virus|viral|phage|coronavirus|influenza|retrovirus|herpesvirus|"
                    r"adenovirus|papillomavirus|norovirus|rotavirus|rhinovirus|"
                    r"paramyxovirus|hepatovirus|poxvirus|flavivirus|filovirus|arenavirus|virus|viridae|virinae|phage|corona|sars|mers|orthomyxo|influenza|adeno|retro|"
                    r"herpes|papilloma|noro|rota|rhino|paramyxo|hepat|pox|flavi|filo|arena|reo|bunya|"
                    r"hiv|ebola|zika|dengue)",
    "fungi": r"(?:fungi|fungal|candida|aspergillus|neurospora|penicillium|trichophyton|"
                    r"cryptococcus|histoplasma|blastomyces|coccidioides|saccharomyces|"
                    r"schizosaccharomyces|pichia|kluyveromyces|mucor|zygomycetes|yeast|myces|mycota|mycet|aspergill|candida|neurospor|penicill|trichophyt|"
                    r"cryptococc|histoplas|blastomy|coccidio|saccharo|schizo|pichia|kluyvero|mucor|zygo)",
    "plants": r"(?:plant|plantae|viridiplantae|rice|wheat|maize|soybean|cotton|populus|"
                    r"arabidopsis|oryza|triticum|zea|glycine|hordeum|barley|avena|oat|"
                    r"sorghum|solanum|potato|tomato|gossypium|brassica|spinacia|"
                    r"spinach|lactuca|lettuce|cucumis|cucumber|melon|vitis|grape|musa|"
                    r"banana|malus|apple|pyrus|pear|prunus|peach|viridiplant|plantae|arabidops|oryza|rice|triticum|wheat|zea|maize|glycine|soybean|"
                    r"hordeum|barley|avena|oat|sorghum|solanum|potato|tomato|gossypium|cotton|populus|"
                    r"brassica|spinacia|spinach|lactuca|lettuce|cucumis|cucumber|melon|vitis|grape|musa|"
                    r"banana|malus|apple|pyrus|pear|prunus|peach)",
    "animals": r"(?:animal|animals|metazoa|metazoan|vertebrate|mammal|mammals|"
                      r"bird|birds|fish|fishes|amphibian|amphibians|reptile|reptiles|insect|insects|"
                      r"arthropod|annelid|mollusk|cnidarian|nematode|echinoderm|porifera|platyhelminth|rotifer|tardigrade|"
                      r"canis|dog|felis|cat|bos|cow|sus|pig|equus|horse|ovis|sheep|capra|goat|oryctolagus|rabbit|macaque|macaca|chimp|"
                      r"pan\s?troglodytes|gorilla|pongo|orangutan|dolphin|delphin|whale|balaen|elephant|loxodonta|monkey)",

    "eukaryotic": r"(?:eukaryotic|eukaryote|eukaryota)",
    "prokaryotic": r"(?:prokaryotic|prokaryote|prokaryota)",
    "synthetic": r"(?:synthetic|vector|plasmid|artificial|construct|chimera|chimeric)",
    "any_rDNA": r"(?:rdna|rrna|ribosomal\s?(dna|rna))"
}


COMPILED = {k: re.compile(v, flags=re.IGNORECASE) for k, v in GROUP_PATTERNS.items()}

SEQTYPE_PATTERNS = {
    "rDNA": r"(?:rdna|rrna|ribosomal\s?(dna|rna))",
    "satellite": r"satellite",
    "repeat": r"repeat",
    "genome": r"(?:whole\s?genome|complete sequence|genome)",
    "telomeric": r"telomer",
    "plasmid": r"plasmid",
    "mitochondrial": r"(?:mitochondria|mitochondrial|mtDNA|mt-?)",
    "chloroplast": r"(?:chloroplast|plastid|cpDNA|cp-[A-Za-z0-9])",
    "contig": r"contig",
    "scaffold": r"scaffold",
    "transposon": r"(?:transposon|retrotransposon|insertion\ssequence|is\d+)",
    "vector": r"(?:cloning\s?vector|expression\s?vector|vector)",
    "synthetic": r"(?:synthetic|artificial|construct|chimera|chimeric)",
    "chromosome": r"(?:chromosome|chr\s?\d+|chromatid|linkage\s?group|lg\s?\d+)",
    "gene": r"(?:gene|cds|coding\s?sequence|open\s?reading\s?frame|orf)",
    "mrna": r"(?:mrna|messenger\s?rna)",
}
SEQTYPE_COMPILED = {k: re.compile(v, re.IGNORECASE) for k, v in SEQTYPE_PATTERNS.items()}

class AutoClustalVisualizer:
    """
    Comprehensive visualizer for AutoClustal results.
    """
    
    def __init__(self, results_dir: str):
        """Initialize visualizer with results directory."""
        self.results_dir = Path(results_dir)
        self.output_dir = self.results_dir / "visualizations"
        self.output_dir.mkdir(exist_ok=True)
        self.organism_categories = {}
        # Model organism categories
        for category in COMPILED:
            self.organism_categories[category] = []  # Initialize list for each category

        
        # {
        #     'Human': ['homo sapiens', 'human', 'h. sapiens', 'homo_sapiens'],
        #     'Mouse': ['mus musculus', 'mouse', 'm. musculus', 'mus_musculus'],
        #     'Rat': ['rattus norvegicus', 'rat', 'r. norvegicus', 'rattus_norvegicus'],
        #     'C. elegans': ['caenorhabditis elegans', 'c. elegans', 'c_elegans', 'nematode'],
        #     'Drosophila': ['drosophila melanogaster', 'd. melanogaster', 'fruit fly', 'drosophila'],
        #     'Zebrafish': ['danio rerio', 'zebrafish', 'd. rerio', 'danio_rerio'],
        #     'Arabidopsis': ['arabidopsis thaliana', 'a. thaliana', 'arabidopsis'],
        #     'Yeast': ['saccharomyces cerevisiae', 's. cerevisiae', 'yeast', 'baker yeast'],
        #     'E. coli': ['escherichia coli', 'e. coli', 'e_coli'],
        #     'Bacteria': ['bacteria', 'bacterial', 'bacillus', 'streptococcus', 'staphylococcus'],
        #     'Virus': ['virus', 'viral', 'phage', 'coronavirus', 'influenza'],
        #     'Fungi': ['fungi', 'fungal', 'candida', 'aspergillus', 'neurospora'],
        #     'Plant': ['plant', 'plantae', 'viridiplantae', 'rice', 'wheat', 'maize'],
        #     'Other': []  # Will be filled with remaining organisms
        # }
    
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

    def _get_seqtype(self, seq_id) -> str:
        s = seq_id.lower()
        for label, pattern in SEQTYPE_COMPILED.items():
            if pattern.search(s):
                return(label)
        return 'unknown'
    
    def _get_canonical(self, seq_id) -> str:
        s = seq_id.lower()
        if re.search(r'\b(human|homo[_\s]?sapiens|h\.?\s?sapiens)\b', s):
            return 'homo sapiens'
        if 'mouse' in s:
            return 'mus musculus'
        if 'rat' in s:
            return 'rattus norvegicus'
        if 'dog' in s:
            return 'canis familiaris'
        if any(x in s for x in ('fish', 'danio', 'rerio', 'zebrafish')):
            return 'danio rerio'
        if 'chicken' in s or 'gallus' in s:
            return 'gallus gallus'
        if 'frog' in s or 'xenopus' in s:
            return 'xenopus laevis'
        if 'yeast' in s or 'saccharomyces' in s:
            return 'saccharomyces cerevisiae'
        if any(x in s for x in ('ecoli', 'coli', 'escherichia')):
            return 'escherichia coli'
        return seq_id

    def _get_latin_name(self, header: str) -> str:
        """
        Extract the Latin binomial name from an NCBI BLAST header or FASTA description.
        Returns 'unknown' if not found.
        """
        # print(f"_get_latin_name: doing {header}")
        if not header:
            # print(f"Header is not present")
            return self._get_canonical(header)

        name = header
        foundmatch = False

        if "no significant" in header.lower():
            # print(f"{header} 0 No significant hits")
            return header
        
        header = re.sub(r'^[A-Z a-z_]+:\s*', '', header)  # Normalize whitespace

        m = re.search(r'^([a-zA-Z]+ [a-zA-Z]+ )', header)
        if m:
            foundmatch = True
            # print(f"{header} 1 Matched species: {m.group(0)}")
            name = m.group(0)

        
        # Case 1: Brackets [Homo sapiens]
        m = re.search(r'\[([A-Z][a-z]+ [a-z]+)\]', header)
        if m:
            foundmatch = True
            # print(f"{header} 1 Matched Genus species: {m.group(1)}")
            name = m.group(1)

        # Case 2: OS=Homo sapiens
        m = re.search(r'OS=([A-Z][a-z]+ [a-z]+)', header)
        if m:
            foundmatch = True
            # print(f"{header} 2 Matched Genus species: {m.group(1)}")
            name = m.group(1)

        m = re.match(r'^([A-Za-z ]+\s+[a-zA-Z]*sp)\.', header)
        if m:
            foundmatch = True
            # print(f"{header} 3 Matched Genus species: {m.group(1)}")
            name = m.group(1)


        # Try Genus species/sp. at beginning of description
        # (stop before "strain|chromosome|plasmid|contig|complete")
        m = re.match(r'([A-Z][a-z]+(?:\s(?:sp\.|[a-z]+)))', header)
        if m:
            foundmatch = True
            # print(f"{header} 3 Matched Genus species: {m.group(1)}")
            name = m.group(1)
    
        name = self._get_canonical(name)
    
        # if not foundmatch:
            # print(f"{header} 4 Can't find!")
            # print(f"_get_canonical: returning {name}")
        
        return(name)

    def _get_category(self, seq_id: str) -> str:
        s = seq_id.lower()
        for category in COMPILED:
            if COMPILED[category].search(s):
                return category
        return "unknown"
    
    def _infer_organism_from_id(self, seq_id: str) -> str:

        """Infer organism from sequence ID using regex patterns."""
        s = seq_id.lower()
        category = self._get_category(s)
        seqtype = self._get_seqtype(s)
        canonical = self._get_latin_name(s)
        return category, canonical, seqtype

    def categorize_organisms(self) -> Dict[str, List[str]]:
        print(self.summary_df.head())
        #                               sequence_id   cluster_id  sequence_length                              description                                   best_match_title  best_match_evalue organism
        # 0  LH00516:376:233737LT3:2:1103:1203:10896  cluster_168              150  LH00516:376:233737LT3:2:1103:1203:10896   Ralstonia sp. CP chromosome 1, complete sequence       1.614140e-68  Unknown
        # 1  LH00516:376:233737LT3:2:1103:3367:10896  cluster_167              150  LH00516:376:233737LT3:2:1103:3367:10896                                No significant hits                NaN  Unknown
        # 2  LH00516:376:233737LT3:2:1103:3515:10896  cluster_184              150  LH00516:376:233737LT3:2:1103:3515:10896  Leifsonia sp. 509MF genome assembly, chromosom...       2.247660e-22  Unknown
        # 3  LH00516:376:233737LT3:2:1103:4864:10896  cluster_183              150  LH00516:376:233737LT3:2:1103:4864:10896  Ramlibacter tataouinensis strain DMF-7 chromos...       1.617000e-30  Unknown
        # 4  LH00516:376:233737LT3:2:1103:5382:10896  cluster_117              150  LH00516:376:233737LT3:2:1103:5382:10896  Eukaryotic synthetic construct chromosome 16 >...       6.863510e-67  Unknown

        """Categorize organisms into model organism categorys and annotate seqtype."""
        organism_counts: Dict[str, List[str]] = {}

        # Ensure required columns exist
        if 'organism' not in self.summary_df.columns:
            self.summary_df['organism'] = 'unknown'
        if 'sequence_id' not in self.summary_df.columns:
            self.summary_df['sequence_id'] = ''
        if 'seqtype' not in self.summary_df.columns:
            self.summary_df['seqtype'] = 'unknown'
        if 'category' not in self.summary_df.columns:
            self.summary_df['category'] = 'Other'

        # Initialize counts container
        for category in COMPILED:
            organism_counts[category] = []
        
        # for category in self.organism_categories:
        #     organism_counts[category] = []
    
        # for category in self.organism_categories:
        #     organism_counts[category] = []
    # 
    #         # print(self.summary_df.head())
        # print(organism_counts.keys())
        # Pass 1: annotate each row with seqtype + category/organism
        updated_rows = []
        for idx, row in self.summary_df.iterrows():
            # seq_id = str(row.get('sequence_id', '') or '')
            seq_id = str(row.get('best_match_title','') or '')
            org_str = str(row.get('organism', 'unknown') or 'unknown').strip()

            # seqtype detection first (independent label)
            # seqtype = self._infer_seqtype(seq_id)
            # organism/category
            category, canonical, seqtype = 'unknown', seq_id, 'unknown'
            if org_str.lower() in ['unknown', 'not searched', 'no significant hits', '', 'nan']:
                category, canonical, seqtype = self._infer_organism_from_id(seq_id)
                # print(f"Inferred {seq_id}: category={category}, canonical={canonical}, seqtype={seqtype})")

            # write back
            self.summary_df.at[idx, 'seqtype'] = seqtype
            # Prefer canonical inferred organism if the original was unknown
            # if org_str.lower() in ['unknown', 'not searched', 'no significant hits', '', 'nan'] and canonical != 'unknown':
            self.summary_df.at[idx, 'organism'] = canonical
            # Always set a category
            self.summary_df.at[idx, 'category'] = category
            if category not in self.organism_categories:
                self.organism_categories[category] = []
            self.organism_categories[category].append(canonical)
            # print(self.summary_df.at[idx, 'organism'])
            # Count by category (store the organism string for your later plots)
            # print(f"here1 appending {category} {canonical}")
            if category not in organism_counts:
                organism_counts[category] = []
            organism_counts[category].append(canonical)
            # .append(self.summary_df.at[idx, 'organism'])
            # print("here")
        
        # print(organism_counts)
        # print(self.organism_categories)
        # print(organism_counts)
        # sys.exit(1)

        # for category in self.summary_df['category'].unique():
        #     organism_counts[category] = 0
        # for row in self.summary_df.iterrows():            
        #     category_str = str(row[1].get('category', '')).strip()
        #     if category_str not in organism_counts:
        #         organism_counts[category_str] = 0
        #     organism_counts[category_str] = organism_counts[category_str] + 1
        # print(organism_counts)
        return organism_counts

    def create_organism_barplot(self):
        """Create comprehensive organism distribution bar plots."""
        print("Creating organism distribution plots...")
        
        # Categorize organisms
        organism_counts = self.categorize_organisms()
        
        # sys.exit(1)
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
    
    def _replace_tree_labels_with_organisms(self):
        """Replace seq_<number> labels in tree with organism names from summary_df."""
        if not hasattr(self, 'summary_df') or self.summary_df is None:
            print("⚠ No summary data available for tree label replacement")
            return self.tree_newick
        
        # Create mapping from sequence index to organism name
        seq_to_organism = {}
        
        # Method 1: Try to extract sequence index from sequence_id if it contains seq_<number>
        for idx, row in self.summary_df.iterrows():
            cluster_id = str(row.get('cluster_id', ''))
            organism = str(row.get('organism', 'unknown'))
            # print(organism)
            # Look for cluster_<number> pattern in cluster_id
            seq_match = re.search(r'cluster_(\d+)', cluster_id)
            if seq_match:
                seq_num = seq_match.group(1)
                seq_to_organism[f"seq_{seq_num}"] = organism
        
        # Method 2: If no seq_<number> found in cluster_id, use dataframe index
        if not seq_to_organism:
            print("⚠ No seq_<number> patterns found in cluster_id, using dataframe index")
            for idx, row in self.summary_df.iterrows():
                organism = str(row.get('organism', 'unknown'))
                seq_to_organism[f"seq_{idx}"] = organism
        
        # Method 3: If organism column has 'unknown' values, try to use best_match_title
        for seq_label, organism in seq_to_organism.items():
            if organism.lower() in ['unknown', 'not searched', 'no significant hits', '', 'nan']:
                # Find the corresponding row and use best_match_title
                seq_num = seq_label.split('_')[1]
                for idx, row in self.summary_df.iterrows():
                    if (f"seq_{idx}" == seq_label or 
                        re.search(rf'seq_{seq_num}', str(row.get('cluster_id', '')))):
                        best_match = str(row.get('best_match_title', 'unknown'))
                        if best_match and best_match.lower() not in ['unknown', 'no significant hits', 'nan']:
                            # Extract organism name from best_match_title
                            organism_from_blast = self._get_latin_name(best_match)
                            seq_to_organism[seq_label] = organism_from_blast
                        break
        
        # Replace seq_<number> labels in tree with organism names
        modified_tree = self.tree_newick
        for seq_label, organism in seq_to_organism.items():
            # Clean up organism name for tree display
            clean_organism = organism.replace(' ', '_').replace('(', '').replace(')', '')
            clean_organism = re.sub(r'[^\w\-_.]', '', clean_organism)  # Remove special chars
            modified_tree = modified_tree.replace(seq_label, clean_organism)
        
        print(f"✓ Replaced {len(seq_to_organism)} sequence labels with organism names")
        return modified_tree

    def create_tree_visualization(self):
        """Create phylogenetic tree visualization."""
        if not self.tree_newick:
            print("⚠ No tree data available")
            return
        
        print("Creating phylogenetic tree visualization...")
        
        # Replace seq_<number> with organism names from summary_df
        tree_with_organisms = self._replace_tree_labels_with_organisms()
        
        try:
            # Try to use BioPython for tree visualization
            from Bio import Phylo
            from io import StringIO
            
            # Parse tree
            tree_handle = StringIO(tree_with_organisms)
            tree = Phylo.read(tree_handle, 'newick')
            
            # Create tree plot
            # print(len(tree.get_terminals()))
            fig, ax = plt.subplots(1, 1, figsize=(24, min(50, len(tree.get_terminals()) / 5)))
            
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
        # Replace seq_<number> with organism names for display
        tree_with_organisms = self._replace_tree_labels_with_organisms()
        
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        
        ax.text(0.5, 0.6, 'Phylogenetic Tree', ha='center', va='center', 
                fontsize=16, fontweight='bold')
        ax.text(0.5, 0.5, f'Newick Format:', ha='center', va='center', 
                fontsize=12, fontweight='bold')
        
        # Display truncated Newick string with organism names
        newick_display = tree_with_organisms[:200] + "..." if len(tree_with_organisms) > 200 else tree_with_organisms
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
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(90, 12))
        fig.suptitle('Detailed Cluster Analysis with Annotations', fontsize=16, fontweight='bold')
        
        # 1. Cluster size distribution with best matches
        labels = []
        short_matches = []
        for cluster_id in cluster_counts.index:
            best_match = cluster_info[cluster_id]['best_match']
            # Custom label extraction logic
            if best_match.lower() == "no significant hits":
                short_match = "Nomatch"
            elif "predicted" in best_match.lower():
                words = best_match.split()
                idx = next((i for i, w in enumerate(words) if w.lower() == "predicted"), None)
                if idx is not None:
                    short_match = " ".join(words[idx:idx+3])
                else:
                    short_match = " ".join(words[:2])
            elif "rDNA" in best_match or "rRNA" in best_match or "rdna" in best_match or "rrna" in best_match:
                short_match = "rRNA"
            elif "homo sapiens" in best_match.lower() or "human" in best_match.lower():
                short_match = "homo sapiens"
            elif ";" in best_match:
                short_match = best_match.split(";")[0]
            elif "MAG:" in best_match:
                mag_match = best_match.split()
                if mag_match:
                    short_match = " ".join(mag_match[:2])
                else:
                    short_match = " ".join(best_match.split()[:2])
            else:
                short_match = " ".join(best_match.split()[:2])
                labels.append(f"{short_match}")
                short_matches.append(short_match)

        # Sort by x axis (short_match) and cluster size
        sorted_items = sorted(zip(labels, cluster_counts.index, cluster_counts.values), key=lambda x: (x[0], -x[2]))
        sorted_labels, sorted_cluster_ids, sorted_counts = zip(*sorted_items) if sorted_items else ([], [], [])

        bars1 = ax1.bar(range(len(sorted_counts)), sorted_counts, 
                    color='skyblue', alpha=0.8, edgecolor='black')
        ax1.set_title('Sequences per Cluster (with Best Match)', fontweight='bold')
        ax1.set_xlabel('Cluster ID + Annotation')
        ax1.set_ylabel('Number of Sequences')
        ax1.set_xticks(range(len(sorted_counts)))
        ax1.set_xticklabels(sorted_labels, rotation=45, ha='right', fontsize=8)

        # Add value labels
        for i, (bar, count) in enumerate(zip(bars1, sorted_counts)):
            ax1.text(bar.get_x() + bar.get_width()/2., count + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')

        # Save sorted plot with width == number of x axis items
        sorted_fig, sorted_ax = plt.subplots(figsize=(max(8, len(sorted_labels)), 6))
        sorted_bars = sorted_ax.bar(range(len(sorted_counts)), sorted_counts, 
                        color='skyblue', alpha=0.8, edgecolor='black')
        sorted_ax.set_title('Sequences per Cluster (Sorted by Annotation)', fontweight='bold')
        sorted_ax.set_xlabel('Cluster ID + Annotation')
        sorted_ax.set_ylabel('Number of Sequences')
        sorted_ax.set_xticks(range(len(sorted_counts)))
        sorted_ax.set_xticklabels(sorted_labels, rotation=90, ha='right', fontsize=8)
        for i, (bar, count) in enumerate(zip(sorted_bars, sorted_counts)):
            sorted_ax.text(bar.get_x() + bar.get_width()/2., count + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')
        plt.tight_layout()
        sorted_path = self.output_dir / 'cluster_analysis_sorted.png'
        sorted_fig.savefig(sorted_path, dpi=300, bbox_inches='tight')
        plt.close(sorted_fig)
        print(f"✓ Sorted cluster analysis saved: {sorted_path}")

        # Group by short_match histogram
        short_match_hist = Counter(short_matches)
        hist_labels, hist_counts = zip(*short_match_hist.items()) if short_match_hist else ([], [])
        hist_fig, hist_ax = plt.subplots(figsize=(max(8, len(hist_labels)), 6))
        hist_bars = hist_ax.bar(range(len(hist_counts)), hist_counts, 
                        color='mediumseagreen', alpha=0.8, edgecolor='black')
        hist_ax.set_title('Clusters Grouped by Annotation (short_match)', fontweight='bold')
        hist_ax.set_xlabel('Annotation Category')
        hist_ax.set_ylabel('Number of Clusters')
        hist_ax.set_xticks(range(len(hist_labels)))
        hist_ax.set_xticklabels(hist_labels, rotation=90, ha='right', fontsize=8)
        for i, (bar, count) in enumerate(zip(hist_bars, hist_counts)):
            hist_ax.text(bar.get_x() + bar.get_width()/2., count + 0.1,
                    str(count), ha='center', va='bottom', fontweight='bold')
        plt.tight_layout()
        hist_path = self.output_dir / 'cluster_annotation_histogram.png'
        hist_fig.savefig(hist_path, dpi=300, bbox_inches='tight')
        plt.close(hist_fig)
        print(f"✓ Cluster annotation histogram saved: {hist_path}")

        # 2. Cluster size histogram
        ax2.hist(cluster_counts.values, bins=max(1, len(cluster_counts)//2), 
                alpha=0.7, color='lightcoral', edgecolor='black')
        ax2.set_title('Cluster Size Distribution', fontweight='bold')
        ax2.set_xlabel('Cluster Size (sequences)')
        ax2.set_ylabel('Number of Clusters')
        
        # 3a. Pie chart of clusters by annotation category (short_match)
        short_match_counts = {}
        for short_match in labels:
            short_match_counts[short_match] = short_match_counts.get(short_match, 0) + 1

        # Only show top N categories, category rest as 'Other'
        N = 8
        sorted_short_matches = sorted(short_match_counts.items(), key=lambda x: x[1], reverse=True)
        top_short_matches = sorted_short_matches[:N]
        other_count = sum([count for _, count in sorted_short_matches[N:]])
        pie_labels = [sm for sm, _ in top_short_matches]
        pie_sizes = [count for _, count in top_short_matches]
        if other_count > 0:
            pie_labels.append('Other')
            pie_sizes.append(other_count)

        colors_pie = plt.cm.Set3(np.linspace(0, 1, len(pie_labels)))
        wedges, texts, autotexts = ax3.pie(pie_sizes, labels=pie_labels, autopct='%1.1f%%',
                    colors=colors_pie, startangle=90)
        ax3.set_title('Clusters by Annotation Category', fontweight='bold')
        # Make percentage text bold
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontweight('bold')

        # Save pie chart separately
        pie_fig, pie_ax = plt.subplots(figsize=(8, 8))
        pie_ax.pie(pie_sizes, labels=pie_labels, autopct='%1.1f%%',
                colors=colors_pie, startangle=90)
        pie_ax.set_title('Clusters by Annotation Category', fontweight='bold')
        pie_path = self.output_dir / 'cluster_annotation_pie.png'
        plt.savefig(pie_path, dpi=300, bbox_inches='tight')
        plt.close(pie_fig)
        print(f"✓ Cluster annotation pie chart saved: {pie_path}")

        # # 3. Annotation success rate per cluster
        # annotated_clusters = 0
        # unannotated_clusters = 0
        
        # for cluster_id, info in cluster_info.items():
        #     if info['best_match'] != "No significant hits":
        #     annotated_clusters += 1
        #     else:
        #     unannotated_clusters += 1
        
        # if annotated_clusters > 0 or unannotated_clusters > 0:
        #     sizes = [annotated_clusters, unannotated_clusters]
        #     labels_pie = ['Annotated Clusters', 'Unannotated Clusters']
        #     colors = ['lightgreen', 'lightcoral']
            
        #     wedges, texts, autotexts = ax3.pie(sizes, labels=labels_pie, autopct='%1.1f%%',
        #                       colors=colors, startangle=90)
        #     ax3.set_title('Cluster Annotation Success Rate', fontweight='bold')
        #     # Make percentage text bold
        #     for autotext in autotexts:
        #     autotext.set_color('white')
        #     autotext.set_fontweight('bold')
        #     # Make percentage text bold
        #     for autotext in autotexts:
        #         autotext.set_color('white')
        #         autotext.set_fontweight('bold')
        
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
            ['Smallest Cluster', smallest_cluster]
            # ['Annotated Clusters', annotated_clusters],
            # ['Unannotated Clusters', unannotated_clusters]
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
        # if 'cluster_id' in self.summary_df.columns:
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
            # f"• Cluster Analysis: {self.output_dir}/cluster_analysis_detailed.png",
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
        # Replace seq_<number> in tree_newick with organism names from summary_df
        # Build mapping: seq_<number> -> organism
        # Create visualizations
        organism_counts = self.create_organism_barplot()
        self.create_tree_visualization()
        # self.create_cluster_analysis()
        
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