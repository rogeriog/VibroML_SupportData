import pandas as pd
import numpy as np
import re
from pathlib import Path

# ---------------- CONFIGURATION ---------------- #
DFT_FILENAME = "neb_image_details_DFT.csv"
MLIP_FOLDER_PATTERN = "LiF_*"
MLIP_FILE_PATTERN = "results_*.csv"
# ----------------------------------------------- #

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', str(s))]

def get_path_stats(energies, forces):
    if len(energies) == 0:
        return None, None, None
    e, f = np.array(energies), np.array(forces)
    return (np.max(e) - e[0]), (e[-1] - e[0]), np.max(f)

def process_dft_file(filepath):
    print(f"--> Processing DFT file: {filepath.name}")
    try:
        df = pd.read_csv(filepath)
        stats = []
        for group_name, group_df in df.groupby("Group"):
            group_df = group_df.copy()
            group_df['sort_key'] = group_df['Image'].apply(natural_sort_key)
            group_df = group_df.sort_values('sort_key')
            b, rxn, f = get_path_stats(
                group_df['Energy_per_atom_eV'].values,
                group_df['Max_Force_eV_A'].values
            )
            stats.append({
                "Group": group_name, "Model": "DFT",
                "Barrier_eV_atom": b, "Rxn_Energy_eV_atom": rxn, "Max_Force_Path": f
            })
        return stats
    except Exception as e:
        print(f"!! Error reading DFT file: {e}")
        return []

def process_mlip_folder(folder_path):
    folder_stats = []
    group_name = folder_path.name
    
    for csv_file in folder_path.glob(MLIP_FILE_PATTERN):
        model_name = csv_file.stem.replace("results_", "").upper().replace("MACE_", "MACE-")
        try:
            df = pd.read_csv(csv_file)
            e_cols = [c for c in df.columns if "energy" in c.lower() and "atom" in c.lower()]
            f_cols = [c for c in df.columns if "force" in c.lower()]
            
            if not e_cols or not f_cols: continue

            if 'filename' in df.columns:
                df['sort_key'] = df['filename'].apply(natural_sort_key)
                df = df.sort_values('sort_key')
            
            b, rxn, f = get_path_stats(df[e_cols[0]].values, df[f_cols[0]].values)
            folder_stats.append({
                "Group": group_name, "Model": model_name,
                "Barrier_eV_atom": b, "Rxn_Energy_eV_atom": rxn, "Max_Force_Path": f
            })
        except Exception as e:
            print(f"   Error in {csv_file.name}: {e}")
    return folder_stats

def calculate_errors_and_aggregate(df):
    """Calculates per-path errors and then aggregates by Model."""
    print("\n--> Calculating Errors vs DFT Reference...")
    
    # 1. Calculate Differences (Model - DFT)
    for col in ["Diff_Barrier", "Diff_Rxn", "Diff_MaxForce"]:
        df[col] = np.nan

    for group in df["Group"].unique():
        dft_row = df[(df["Group"] == group) & (df["Model"] == "DFT")]
        if dft_row.empty: continue
            
        ref_b = dft_row.iloc[0]["Barrier_eV_atom"]
        ref_r = dft_row.iloc[0]["Rxn_Energy_eV_atom"]
        ref_f = dft_row.iloc[0]["Max_Force_Path"]
        
        mask = df["Group"] == group
        df.loc[mask, "Diff_Barrier"] = df.loc[mask, "Barrier_eV_atom"] - ref_b
        df.loc[mask, "Diff_Rxn"] = df.loc[mask, "Rxn_Energy_eV_atom"] - ref_r
        df.loc[mask, "Diff_MaxForce"] = df.loc[mask, "Max_Force_Path"] - ref_f

    # 2. Aggregation: Calculate Mean Absolute Error (MAE) per Model
    # Filter out DFT rows for the stats summary so we only judge MLIPs
    mlip_df = df[df["Model"] != "DFT"].copy()
    
    # We take the absolute value of differences for MAE
    mlip_df["Abs_Diff_Barrier"] = mlip_df["Diff_Barrier"].abs()
    mlip_df["Abs_Diff_Rxn"] = mlip_df["Diff_Rxn"].abs()
    mlip_df["Abs_Diff_Force"] = mlip_df["Diff_MaxForce"].abs()
    
    agg_stats = mlip_df.groupby("Model").agg({
        "Abs_Diff_Barrier": "mean",
        "Abs_Diff_Rxn": "mean",
        "Abs_Diff_Force": "mean"
    }).reset_index()
    
    agg_stats.columns = ["Model", "MAE_Barrier", "MAE_Rxn_Energy", "MAE_Max_Force"]
    
    return df, agg_stats

def main():
    base_path = Path(".")
    all_data = []

    # 1. Read Data
    dft_file = base_path / DFT_FILENAME
    if dft_file.exists():
        all_data.extend(process_dft_file(dft_file))
    
    for subdir in sorted([d for d in base_path.glob(MLIP_FOLDER_PATTERN) if d.is_dir()]):
        all_data.extend(process_mlip_folder(subdir))

    if not all_data:
        print("No data found.")
        return

    # 2. Process
    final_df = pd.DataFrame(all_data)
    final_df, agg_df = calculate_errors_and_aggregate(final_df)

    # 3. Save & Print Detailed
    cols = ["Group", "Model", "Barrier_eV_atom", "Diff_Barrier", 
            "Rxn_Energy_eV_atom", "Diff_Rxn", "Max_Force_Path", "Diff_MaxForce"]
    final_df = final_df[cols].sort_values(['Group', 'Model'])
    final_df.to_csv("final_path_statistics_detailed.csv", index=False, float_format="%.6f")

    # 4. Save & Print Aggregated
    agg_df.to_csv("mlip_error_overview.csv", index=False, float_format="%.6f")

    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    
    print("-" * 80)
    print("DETAILED RESULTS (saved to final_path_statistics_detailed.csv)")
    print("-" * 80)
    print(final_df)
    
    print("\n" + "=" * 80)
    print("MLIP AGGREGATED ERROR OVERVIEW (MAE) (saved to mlip_error_overview.csv)")
    print("MAE = Mean Absolute Error across all paths (Cm, P1, Pm)")
    print("=" * 80)
    print(agg_df)

if __name__ == "__main__":
    main()