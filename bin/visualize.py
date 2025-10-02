#!/usr/bin/env python3
"""
AutoClustal Quick Visualizer
============================

Quick launcher for creating comprehensive visualizations from AutoClustal results.

Created by Stella R. Hartono 2025

"""

import os
import sys
import argparse
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bin import grapher

def main():
    print("🧬 AutoClustal Visualization Tool")
    print("=" * 40)
    
    # Check if results directory is provided
    if len(sys.argv) < 2:
        print("Usage examples:")
        print("  python visualize.py results/")
        print("  python visualize.py results_blast/")
        print("  python visualize.py demo_results/")
        print()
        
        # Look for results directories
        potential_dirs = []
        for item in os.listdir('.'):
            if os.path.isdir(item) and ('result' in item.lower() or item == 'output'):
                potential_dirs.append(item)
        
        if potential_dirs:
            print("Found these results directories:")
            for i, dir_name in enumerate(potential_dirs, 1):
                print(f"  {i}. {dir_name}")
            print()
            
            try:
                choice = input("Enter number to visualize (or press Enter to cancel): ")
                if choice.strip():
                    selected_dir = potential_dirs[int(choice) - 1]
                    print(f"Selected: {selected_dir}")
                else:
                    sys.exit(0)
            except (ValueError, IndexError):
                print("Invalid selection.")
                sys.exit(1)
        else:
            print("No AutoClustal results directories found.")
            print("Run AutoClustal first: python autoclustal.py -i sequences.fasta -o results/")
            sys.exit(1)
    else:
        selected_dir = sys.argv[1]
    
    # Check if directory exists
    if not os.path.exists(selected_dir):
        print(f"❌ Directory not found: {selected_dir}")
        sys.exit(1)
    
    # Check for required files
    summary_file = Path(selected_dir) / "summary_table.csv"
    if not summary_file.exists():
        print(f"❌ Summary table not found: {summary_file}")
        print("This doesn't appear to be an AutoClustal results directory.")
        sys.exit(1)
    
    print(f"📊 Visualizing results from: {selected_dir}")
    print()
    
    # Import and run the grapher
    try:
        # from autoclustal.bin.grapher import AutoClustalVisualizer
        
        visualizer = grapher.AutoClustalVisualizer(selected_dir)
        visualizer.run_visualization()
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("Make sure grapher.py is in the same directory.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Visualization failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()