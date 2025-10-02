#!/usr/bin/env python3
"""
Create Demo Data for AutoClustal Visualization
==============================================

Creates realistic demo data with organism annotations for testing the grapher.

Created by Stella R. Hartono 2025

"""

import pandas as pd
from pathlib import Path
import random

def create_demo_data(output_dir="demo_results"):
    """Create demo data with realistic organism annotations."""
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create demo sequence data with realistic organism distribution
    demo_data = [
        # Human sequences
        {"sequence_id": "human_gene1", "cluster_id": "cluster_0", "sequence_length": 1245, "description": "Human GAPDH gene", "best_match_title": "glyceraldehyde-3-phosphate dehydrogenase [Homo sapiens]", "organism": "Homo sapiens"},
        {"sequence_id": "human_gene2", "cluster_id": "cluster_0", "sequence_length": 1189, "description": "Human actin gene", "best_match_title": "actin beta [Homo sapiens]", "organism": "Homo sapiens"},
        {"sequence_id": "human_gene3", "cluster_id": "cluster_1", "sequence_length": 987, "description": "Human insulin gene", "best_match_title": "insulin [Homo sapiens]", "organism": "Homo sapiens"},
        
        # Mouse sequences
        {"sequence_id": "mouse_gene1", "cluster_id": "cluster_0", "sequence_length": 1234, "description": "Mouse GAPDH gene", "best_match_title": "glyceraldehyde-3-phosphate dehydrogenase [Mus musculus]", "organism": "Mus musculus"},
        {"sequence_id": "mouse_gene2", "cluster_id": "cluster_0", "sequence_length": 1198, "description": "Mouse actin gene", "best_match_title": "actin beta [Mus musculus]", "organism": "Mus musculus"},
        
        # C. elegans sequences
        {"sequence_id": "celegans_gene1", "cluster_id": "cluster_2", "sequence_length": 876, "description": "C. elegans unc gene", "best_match_title": "uncoordinated protein [Caenorhabditis elegans]", "organism": "Caenorhabditis elegans"},
        {"sequence_id": "celegans_gene2", "cluster_id": "cluster_2", "sequence_length": 923, "description": "C. elegans dpy gene", "best_match_title": "dumpy protein [Caenorhabditis elegans]", "organism": "Caenorhabditis elegans"},
        
        # Drosophila sequences
        {"sequence_id": "drosophila_gene1", "cluster_id": "cluster_3", "sequence_length": 1156, "description": "Drosophila white gene", "best_match_title": "white protein [Drosophila melanogaster]", "organism": "Drosophila melanogaster"},
        
        # Yeast sequences
        {"sequence_id": "yeast_gene1", "cluster_id": "cluster_4", "sequence_length": 789, "description": "Yeast ADH gene", "best_match_title": "alcohol dehydrogenase [Saccharomyces cerevisiae]", "organism": "Saccharomyces cerevisiae"},
        {"sequence_id": "yeast_gene2", "cluster_id": "cluster_4", "sequence_length": 834, "description": "Yeast HIS gene", "best_match_title": "histidine biosynthesis [Saccharomyces cerevisiae]", "organism": "Saccharomyces cerevisiae"},
        {"sequence_id": "yeast_gene3", "cluster_id": "cluster_4", "sequence_length": 756, "description": "Yeast TRP gene", "best_match_title": "tryptophan biosynthesis [Saccharomyces cerevisiae]", "organism": "Saccharomyces cerevisiae"},
        
        # E. coli sequences
        {"sequence_id": "ecoli_gene1", "cluster_id": "cluster_5", "sequence_length": 1023, "description": "E. coli lacZ gene", "best_match_title": "beta-galactosidase [Escherichia coli]", "organism": "Escherichia coli"},
        {"sequence_id": "ecoli_gene2", "cluster_id": "cluster_5", "sequence_length": 945, "description": "E. coli ampR gene", "best_match_title": "ampicillin resistance [Escherichia coli]", "organism": "Escherichia coli"},
        
        # Virus sequences
        {"sequence_id": "virus_gene1", "cluster_id": "cluster_6", "sequence_length": 1789, "description": "HIV gag gene", "best_match_title": "gag polyprotein [Human immunodeficiency virus 1]", "organism": "Human immunodeficiency virus 1"},
        {"sequence_id": "virus_gene2", "cluster_id": "cluster_6", "sequence_length": 1234, "description": "Influenza HA gene", "best_match_title": "hemagglutinin [Influenza A virus]", "organism": "Influenza A virus"},
        
        # Plant sequences
        {"sequence_id": "arabidopsis_gene1", "cluster_id": "cluster_7", "sequence_length": 1345, "description": "Arabidopsis RuBisCO gene", "best_match_title": "ribulose bisphosphate carboxylase [Arabidopsis thaliana]", "organism": "Arabidopsis thaliana"},
        {"sequence_id": "rice_gene1", "cluster_id": "cluster_7", "sequence_length": 1289, "description": "Rice amylase gene", "best_match_title": "alpha-amylase [Oryza sativa]", "organism": "Oryza sativa"},
        
        # Zebrafish sequences
        {"sequence_id": "zebrafish_gene1", "cluster_id": "cluster_8", "sequence_length": 1456, "description": "Zebrafish pax gene", "best_match_title": "paired box protein [Danio rerio]", "organism": "Danio rerio"},
        
        # Other organisms
        {"sequence_id": "fungi_gene1", "cluster_id": "cluster_9", "sequence_length": 967, "description": "Aspergillus gene", "best_match_title": "hypothetical protein [Aspergillus niger]", "organism": "Aspergillus niger"},
        {"sequence_id": "bacteria_gene1", "cluster_id": "cluster_10", "sequence_length": 834, "description": "Bacillus gene", "best_match_title": "sporulation protein [Bacillus subtilis]", "organism": "Bacillus subtilis"},
    ]
    
    # Create DataFrame and save as CSV
    df = pd.DataFrame(demo_data)
    csv_path = output_path / "summary_table.csv"
    df.to_csv(csv_path, index=False)
    
    # Create a simple tree file
    tree_newick = "((((human_gene1:0.1,human_gene2:0.1):0.2,(mouse_gene1:0.15,mouse_gene2:0.15):0.2):0.3,((celegans_gene1:0.2,celegans_gene2:0.2):0.25,drosophila_gene1:0.45):0.3):0.4,(((yeast_gene1:0.1,yeast_gene2:0.1,yeast_gene3:0.1):0.3,(ecoli_gene1:0.2,ecoli_gene2:0.2):0.3):0.4,((virus_gene1:0.3,virus_gene2:0.3):0.2,(arabidopsis_gene1:0.25,rice_gene1:0.25):0.2):0.4):0.5);"
    
    tree_path = output_path / "tree.newick"
    with open(tree_path, 'w') as f:
        f.write(tree_newick)
    
    print(f"✓ Demo data created in: {output_path}")
    print(f"  - {len(demo_data)} sequences across {len(df['cluster_id'].unique())} clusters")
    print(f"  - Organisms: {len(df['organism'].unique())} species")
    print(f"✓ Files created:")
    print(f"  - {csv_path}")
    print(f"  - {tree_path}")
    
    return output_path

if __name__ == "__main__":
    demo_dir = create_demo_data()
    
    print(f"\nTo visualize the demo data, run:")
    print(f"python grapher.py {demo_dir}")