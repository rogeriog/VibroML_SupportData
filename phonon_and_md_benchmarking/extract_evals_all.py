import os
import json
import numpy as np
import pandas as pd
from ase.io import read
from ase.formula import Formula
import re
import math

# --- Configuration ---
ENGINE_DIRS = {
    'MACE_evals': 'MACE',
    'MACEMP_evals': 'MACEMP',
    'ESEN_evals': 'ESEN',
    'UMA_evals': 'UMA'
}

ENGINE_SORT_ORDER = ['MACE', 'MACEMP', 'ESEN', 'UMA']

# --- Standard Phase Reference Map ---
# Maps reduced compositions to the specific Reference Run ID.
STANDARD_PHASE_MAP = {
    'HfO2': 'HfO2_Fm3m',
    'FLi': 'LiFsimplecubic',      # ASE reduces LiF -> FLi
    'CsI3Pb': 'CsPbI3',           # ASE reduces CsPbI3 -> CsI3Pb
    'Bi2O7Sn2': 'Sn2Bi2O7'        # ASE reduces Sn2Bi2O7 -> Bi2O7Sn2
}

# --- Parsing Helpers ---

def get_reduced_formula(formula_str):
    """Reduces and sorts formula (e.g. Sn16Bi16O56 -> Bi2O7Sn2)."""
    try:
        f = Formula(formula_str)
        count_dict = f.count()
        if not count_dict: return formula_str
        counts = list(count_dict.values())
        common_divisor = math.gcd(*counts)
        reduced_parts = []
        for sym in sorted(count_dict.keys()):
            count = count_dict[sym] // common_divisor
            reduced_parts.append(sym if count == 1 else f"{sym}{count}")
        return "".join(reduced_parts)
    except:
        return formula_str

def parse_folder_name_info(folder_name):
    """Parses folder name for MACE Energy."""
    info = {'MACE Energy (eV/atom)': None}
    energy_match = re.search(r"_energy_([mp]\d+p\d+)", folder_name)
    if not energy_match:
        energy_match = re.search(r"_([mp]\d+p\d+)(?:_UMA|_ESEN|_MACE|_phonon|$)", folder_name)

    if energy_match:
        try:
            val = energy_match.group(1).replace('m', '-').replace('p', '.')
            info['MACE Energy (eV/atom)'] = float(val)
        except ValueError:
            pass
    return info

def clean_run_id(folder_name):
    """Extracts a normalized Run ID."""
    clean = re.sub(r"^EVAL_cifs_vibroml\d*_[a-zA-Z]+_", "", folder_name)
    clean = re.sub(r"_phonon_output_.*$", "", clean)
    
    for tag in ['_MACE', '_MACEMP', '_ESEN', '_UMA']:
        if clean.endswith(tag):
            clean = clean[:-len(tag)]

    energy_match = re.search(r"_m?\d+p\d+", clean)
    if energy_match:
        clean = clean.split(energy_match.group(0))[0]

    return clean.strip('_')

def get_relaxed_info(run_dir):
    """Reads relaxed structure info."""
    relax_dir = os.path.join(run_dir, "initial_relaxation_for_single_run")
    info = {
        'Composition': "N/A", 'Cell Params': "N/A", 
        'Relaxed Energy (eV/atom)': None, 'Initial Space Group': "N/A", 
        'Phase Retained': "Unknown"
    }
    
    cif_files = [f for f in os.listdir(relax_dir) if f.endswith("_relaxed.cif")] if os.path.exists(relax_dir) else []
    if cif_files:
        try:
            atoms = read(os.path.join(relax_dir, cif_files[0]))
            info['Composition'] = get_reduced_formula(atoms.get_chemical_formula())
            cp = atoms.get_cell().cellpar()
            info['Cell Params'] = f"{cp[0]:.2f}, {cp[1]:.2f}, {cp[2]:.2f}, {cp[3]:.1f}, {cp[4]:.1f}, {cp[5]:.1f}"
        except: pass

    e_file = os.path.join(relax_dir, "energy_info.txt")
    if os.path.exists(e_file):
        with open(e_file, 'r') as f:
            for line in f:
                if "Final energy per atom:" in line:
                    try: info['Relaxed Energy (eV/atom)'] = float(line.split(":")[1].strip().split()[0])
                    except: pass

    def get_sym(fname):
        p = os.path.join(relax_dir, fname)
        if os.path.exists(p):
            with open(p, 'r') as f:
                m = re.search(r"International symbol:\s+(\S+)", f.read())
                if m: return m.group(1)
        return "N/A"

    init_sg = get_sym("initial_symmetry_analysis.txt")
    final_sg = get_sym("relaxed_symmetry_analysis.txt")
    info['Initial Space Group'] = init_sg
    if init_sg != "N/A" and final_sg != "N/A":
        info['Phase Retained'] = "Yes" if init_sg == final_sg else f"No ({final_sg})"

    return info

def get_phonon_info(run_dir):
    """Extracts phonon stability data."""
    info = {
        'Softest Freq': "N/A", 'Highest Freq': "N/A", 
        'Stable': "N/A", 'Softness Fraction': "N/A",
        'Soft HS Point': "N/A"
    }
    
    bs_files = [f for f in os.listdir(run_dir) if f.startswith("band_structure_energies_")]
    if bs_files:
        try:
            bs = np.loadtxt(os.path.join(run_dir, bs_files[0]))
            min_f = np.min(bs)
            all_freqs = bs.flatten()
            total_modes = all_freqs.size
            unstable_count = np.sum(all_freqs < -0.05)
            fraction = unstable_count / total_modes if total_modes > 0 else 0.0
            
            info['Softest Freq'] = round(min_f, 4)
            info['Highest Freq'] = round(np.max(bs), 4)
            info['Softness Fraction'] = f"{fraction:.4f}"

            freq_ok = min_f > -1.5
            frac_ok = fraction < 0.02
            perfect = min_f > -0.05

            if perfect:
                info['Stable'] = "✅ Stable"
            elif freq_ok and frac_ok:
                info['Stable'] = "✅ Stable"
            elif freq_ok or frac_ok:
                info['Stable'] = "⚠️ Quasi-Stable"
            else:
                info['Stable'] = "❌ Unstable"
        except: pass

    spa_path = os.path.join(run_dir, "special_point_analysis.json")
    if os.path.exists(spa_path):
        try:
            with open(spa_path, 'r') as f:
                data = json.load(f)
                softest = min(data, key=lambda x: x['min_frequency'])
                if softest['min_frequency'] < -0.05:
                    info['Soft HS Point'] = softest['label']
                else:
                    info['Soft HS Point'] = "None"
        except: pass
    return info

def process_run(path, engine):
    """Extracts all data for a single run."""
    folder_name = os.path.basename(path.rstrip('/'))
    data = {'Engine': engine, 'Run ID': clean_run_id(folder_name)}
    
    data.update(parse_folder_name_info(folder_name))
    data.update(get_relaxed_info(path))
    
    if data['MACE Energy (eV/atom)'] is None:
        data['MACE Energy (eV/atom)'] = data['Relaxed Energy (eV/atom)']

    if not data['Run ID'] or data['Run ID'] == data.get('Composition', ''):
        comp = data.get('Composition', 'Unknown')
        data['Run ID'] = f"{comp}_Eval"

    data.update(get_phonon_info(path))
    return data

def report_missing_runs(df, expected_engines=None):
    if expected_engines is None:
        expected_engines = ENGINE_SORT_ORDER

    print(f"\n--- Checking for Missing Runs (Expected: {', '.join(expected_engines)}) ---")
    grouped = df.groupby(['Run ID', 'Composition'])['Engine'].apply(set).reset_index()
    
    missing_data = []
    expected_set = set(expected_engines)

    for index, row in grouped.iterrows():
        present = row['Engine']
        missing = expected_set - present
        if missing:
            missing_data.append({
                'Run ID': row['Run ID'],
                'Composition': row['Composition'],
                'Missing Engines': ", ".join(sorted(list(missing))),
                'Missing Count': len(missing)
            })

    if not missing_data:
        print("✅ Complete! All Run IDs have results for all engines.")
        return

    df_missing = pd.DataFrame(missing_data)
    df_missing = df_missing.sort_values(by=['Composition', 'Run ID'])
    
    print(f"⚠️  Found {len(df_missing)} incomplete samples.")
    print(df_missing.to_markdown(index=False))
    
    csv_out = "vibroml_missing_runs.csv"
    df_missing.to_csv(csv_out, index=False)
    print(f"\nMissing runs report saved to: {csv_out}")

# --- Main Execution ---

all_data = []
cwd = os.getcwd()

print(f"Scanning directories in {cwd}...")

# 1. Collect Data
for dir_name, engine_name in ENGINE_DIRS.items():
    full_dir_path = os.path.join(cwd, dir_name)
    if os.path.exists(full_dir_path):
        print(f"Processing {engine_name} runs in {dir_name}...")
        for run_folder in os.listdir(full_dir_path):
            run_path = os.path.join(full_dir_path, run_folder)
            if os.path.isdir(run_path) and "phonon_output" in run_folder:
                try:
                    all_data.append(process_run(run_path, engine_name))
                except Exception as e:
                    print(f"  Error reading {run_folder}: {e}")

if not all_data:
    print("No data found. Exiting.")
    exit()

# 2. Create DataFrame
df = pd.DataFrame(all_data)

# Deduplication
df['is_na'] = df['Stable'].apply(lambda x: 1 if x == "N/A" else 0)
df = df.sort_values(by=['Engine', 'Run ID', 'is_na'], ascending=[True, True, True])
df = df.drop_duplicates(subset=['Engine', 'Run ID'], keep='first').drop(columns=['is_na'])

# --- NEW: Calculate Relative Energy (meV/atom) ---

def calculate_relative_energies(df, ref_map):
    """
    Calculates Relative Relaxed Energy (meV/atom) based on standard references.
    Formula: (Sample Energy - Reference Energy) * 1000
    """
    # 1. Build Reference Lookup Dictionary: (Engine, Composition) -> Energy (eV)
    ref_energies = {}
    
    for engine in df['Engine'].unique():
        for comp, ref_run_id in ref_map.items():
            # Find the exact row matching Engine, Composition, and the Reference Run ID
            ref_row = df[
                (df['Engine'] == engine) & 
                (df['Composition'] == comp) & 
                (df['Run ID'] == ref_run_id)
            ]
            
            if not ref_row.empty:
                val = ref_row['Relaxed Energy (eV/atom)'].values[0]
                if val is not None and not np.isnan(val):
                    ref_energies[(engine, comp)] = val
    
    # 2. Apply calculation and conversion
    def get_relative(row):
        key = (row['Engine'], row['Composition'])
        val = row['Relaxed Energy (eV/atom)']
        
        # Check if we have a reference energy for this composition/engine combo
        if key in ref_energies and val is not None:
            diff_ev = val - ref_energies[key]
            # Convert eV to meV (multiply by 1000) and round to 1 decimal
            return round(diff_ev * 1000.0, 1)
        return None

    # New column name reflects unit
    col_name = 'Relative Relaxed Energy (meV/atom)'
    df[col_name] = df.apply(get_relative, axis=1)
    return df

print("\nCalculating relative energies (meV/atom)...")
df = calculate_relative_energies(df, STANDARD_PHASE_MAP)

# -----------------------------------------------------

# 3. Handle Sorting Energy
mace_rows = df[df['Engine'] == 'MACE']
id_to_energy = mace_rows.set_index('Run ID')['MACE Energy (eV/atom)'].to_dict()

def get_sort_energy(row):
    val = id_to_energy.get(row['Run ID'])
    if val is None: val = row['MACE Energy (eV/atom)']
    return val if val is not None else -9999.0

df['Sort_Energy'] = df.apply(get_sort_energy, axis=1)
df['Engine'] = pd.Categorical(df['Engine'], categories=ENGINE_SORT_ORDER, ordered=True)

# Final Sort for Display
df = df.sort_values(
    by=['Composition', 'Sort_Energy', 'Run ID', 'Engine'], 
    ascending=[True, False, True, True]
)

# 4. Final Columns 
final_cols = [
    'Engine', 'Run ID', 'Composition', 'Initial Space Group', 
    'Relative Relaxed Energy (meV/atom)',  # Updated column name
    'Stable', 'Softness Fraction', 
    'Softest Freq', 'Soft HS Point', 
    'Relaxed Energy (eV/atom)', 
    'Phase Retained', 'Cell Params'
]
df_final = df[[c for c in final_cols if c in df.columns]]

# 5. Outputs
print("\n--- Final Consolidated Summary (First 20 Rows) ---")
print(df_final.head(20).to_markdown(index=False))

csv_name = "vibroml_master_eval_summary.csv"
df_final.to_csv(csv_name, index=False)
print(f"\nFull results saved to: {csv_name}")

# 6. Run Missing Check
report_missing_runs(df)