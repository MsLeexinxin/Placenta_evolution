#!/usr/bin/env python3

import pandas as pd
import os
import sys
import json

def main():
    if len(sys.argv) != 5:
        print("Error: Missing arguments.")
        print("Usage: python 05_filter_cpdb_results.py <input_dir> <cellpairs.txt> <groups.json> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    cellpairs_file = sys.argv[2]
    groups_file = sys.argv[3]
    output_dir = sys.argv[4]

    os.makedirs(output_dir, exist_ok=True)

    
    # Load grouping json (e.g., {"Hemo": ["H", "M", "MA"], "Endo": ["SQ", "D"]})
    try:
        with open(groups_file, 'r') as f:
            groups = json.load(f)
        print(f"Loaded {len(groups)} groups from {groups_file}")
    except Exception as e:
        print(f"Error loading {groups_file}: {e}")
        sys.exit(1)

    # Load target cell-cell pairs
    try:
        with open(cellpairs_file, 'r') as f:
            cellpairs = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(cellpairs)} cell pairs to analyze.")
    except Exception as e:
        print(f"Error loading {cellpairs_file}: {e}")
        sys.exit(1)

    # Loop through cell pairs
    for cp in cellpairs:
        interaction_metadata = {}
        
        # Collect all interaction metadata of the current cell pair
        for group, species_list in groups.items():
            for sp in species_list:
                file_path = os.path.join(input_dir, f"{sp}-pvalues.txt")
                if not os.path.exists(file_path):
                    continue
                try:
                    df = pd.read_csv(file_path, sep="\t")
                    if cp not in df.columns:
                        continue
                    for _, row in df.iterrows():
                        interaction = row["interacting_pair"]
                        if interaction not in interaction_metadata:
                            interaction_metadata[interaction] = {
                                "gene_a": row["gene_a"],
                                "gene_b": row["gene_b"]
                            }
                except Exception as e:
                    print(f"Warning: Could not process {file_path} for metadata: {e}")
                    continue

        if not interaction_metadata:
            continue

        # Initialize data(default p=1.0)
        all_data = {}
        for interaction, meta in interaction_metadata.items():
            all_data[interaction] = {
                "gene_a": meta["gene_a"],
                "gene_b": meta["gene_b"]
            }
            for group, species_list in groups.items():
                all_data[interaction][group] = {sp: 1.0 for sp in species_list}

        # Extract p-values
        for group, species_list in groups.items():
            for sp in species_list:
                file_path = os.path.join(input_dir, f"{sp}-pvalues.txt")
                if not os.path.exists(file_path):
                    continue
                try:
                    df = pd.read_csv(file_path, sep="\t")
                    if cp not in df.columns:
                        continue
                    
                    interaction_pvals = dict(zip(df["interacting_pair"], df[cp]))
                    
                    # p < 0.05 for all species within the group
                    for interaction in all_data:
                        if interaction in interaction_pvals:
                            all_data[interaction][group][sp] = float(interaction_pvals[interaction])
                except Exception as e:
                    continue

        # Validate the conditions & save the results 
        for group in groups:
            output = []
            for interaction, data in all_data.items():
                current_p = [data[group][sp] for sp in groups[group]]
                
                if all(p < 0.05 for p in current_p):
                    output.append({
                        "interacting_pair": interaction,
                        "gene_a": data["gene_a"],
                        "gene_b": data["gene_b"]
                    })
            
            if output:
                df_out = pd.DataFrame(output)
                safe_cp = cp.replace("|", "_")
                out_file = os.path.join(output_dir, f"{group}-{safe_cp}.result.csv")
                df_out.to_csv(out_file, index=False)
                print(f"  -> Saved significant interactions for {group} in {safe_cp}")


if __name__ == "__main__":
    main()