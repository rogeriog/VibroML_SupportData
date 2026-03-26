import os
import re
import json
import numpy as np
import pandas as pd
from ase.io import read
from ase.formula import Formula
import math
import glob

# --- Configuration ---
ENGINE_DIRS = {
    'MACE_MD_evals': 'MACE OMAT',
    'MACEMP_MD_evals': 'MACE MP',
    'ESEN_MD_evals': 'eSEN',
    'UMA_MD_evals': 'UMA'
}

ENGINE_SORT_ORDER = ['MACE OMAT', 'MACE MP', 'eSEN', 'UMA']
COMPOUND_ORDER = ['LiF', 'CsPbI3', 'HfO2', 'Bi2Sn2O7']

# CIF source directory for comparison
ROOT_CIF_DIR = "/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/cifs_to_eval"

# --- Parsing Helpers ---

def get_reduced_formula(formula_str):
    # NOTE: This function produces alphabetically sorted formulas (e.g., LiF -> FLi).
    # We will normalize this back to the user's preferred order in main().
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

def parse_energy_from_folder(folder_name):
    match = re.search(r'(?:_|energy_)([mp])(\d+)p(\d+)', folder_name)
    if match:
        sign_char = match.group(1)
        integer_part = match.group(2)
        decimal_part = match.group(3)
        value = float(f"{integer_part}.{decimal_part}")
        if sign_char == 'm': value = -value
        return value
    return None

def parse_folder_name_info(folder_name):
    info = {'Energy (eV/atom)': None}
    info['Energy (eV/atom)'] = parse_energy_from_folder(folder_name)
    return info

def clean_run_id(folder_name):
    clean = re.sub(r'^MD_cifs_vibroml\d*_(?:mace|macemp|esen|uma)_', '', folder_name.lower())
    clean = re.sub(r'^MD_cifs_', '', clean)
    clean = re.sub(r'_(?:mace|macemp|esen|uma)_md_stability.*$', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'_(?:mace|macemp|esen|uma)_md.*$', '', clean, flags=re.IGNORECASE)
    clean = re.sub(r'_\d{8}-\d{6}$', '', clean)
    clean = re.sub(r'_phonon_output$', '', clean)
    
    specific_id_match = re.search(
        r'((?:optrandom_|opt_random_)?(?:unique|top|targeted_ga|trad|trad_all|traditional_all_ga)_?\d*(?:_sg_[a-z0-9]+)?(?:_iter\d+)?(?:_sample\d+)?)',
        clean, re.IGNORECASE
    )
    
    if specific_id_match:
        return specific_id_match.group(1)
    return clean

def get_relaxed_info(run_dir):
    relax_dir = os.path.join(run_dir, "initial_relaxation_for_single_run")
    info = {
        'Composition': "N/A", 'Cell Params': "N/A", 
        'Relaxed Energy (eV/atom)': None, 'Initial Space Group': "N/A", 
        'Phase Retained': "Unknown"
    }
    
    # 1. READ CIF for Composition
    if os.path.exists(relax_dir):
        cif_files = [f for f in os.listdir(relax_dir) if f.endswith("_relaxed.cif")]
        if cif_files:
            cif_path = os.path.join(relax_dir, cif_files[0])
            try:
                atoms = read(cif_path)
                # Get formula and reduce it (e.g., Hf4O8 -> HfO2)
                raw_formula = atoms.get_chemical_formula()
                info['Composition'] = get_reduced_formula(raw_formula)
                
                cp = atoms.get_cell().cellpar()
                info['Cell Params'] = f"{cp[0]:.2f}, {cp[1]:.2f}, {cp[2]:.2f}, {cp[3]:.1f}, {cp[4]:.1f}, {cp[5]:.1f}"
            except Exception as e:
                print(f"⚠️ Error reading CIF at {cif_path}: {e}")

    # 2. READ Energy Info
    e_file = os.path.join(relax_dir, "energy_info.txt")
    if os.path.exists(e_file):
        with open(e_file, 'r') as f:
            for line in f:
                if "Final energy per atom:" in line:
                    try: info['Relaxed Energy (eV/atom)'] = float(line.split(":")[1].strip().split()[0])
                    except: pass

    # 3. READ Symmetry Info
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

def parse_md_report(report_path):
    data = {
        "verdict": "N/A", "confidence": "N/A",
        "rmsd_val": "N/A", "rmsd_verdict": "N/A",
        "vol_val": "N/A", "vol_verdict": "N/A",
        "rdf_val": "N/A", "rdf_verdict": "N/A",
        "sym_verdict": "N/A", "sym_group": "N/A", "post_md_sgs": "N/A",
        "status": "COMPLETED"
    }
    
    if not os.path.exists(report_path):
        data["status"] = "MISSING_REPORT"
        return data

    try:
        with open(report_path, 'r') as f: 
            content = f.read()
            lines = content.splitlines()
        
        if "FAILED" in content or "Error" in content or "error" in content:
            data["status"] = "FAILED"
            
        v_match = re.search(r"Stability Verdict:\s+(\w+)", content)
        c_match = re.search(r"Confidence Level:\s+(\w+)", content)
        data["verdict"] = v_match.group(1) if v_match else "N/A"
        data["confidence"] = c_match.group(1) if c_match else "N/A"

        # --- Metrics Parsing ---
        m_rmsd = re.search(r"RMSD:\s+([\d\.]+)\s+Å.*-\s+(PASS|FAIL)", content)
        if m_rmsd:
            data['rmsd_val'] = float(m_rmsd.group(1))
            data['rmsd_verdict'] = m_rmsd.group(2)

        m_vol = re.search(r"Volume change:\s+([\d\.]+)\%.*-\s+(PASS|FAIL)", content)
        if m_vol:
            data['vol_val'] = float(m_vol.group(1))
            data['vol_verdict'] = m_vol.group(2)

        m_rdf = re.search(r"RDF correlation:\s+([-\d\.]+).*-\s+(PASS|FAIL)", content)
        if m_rdf:
            data['rdf_val'] = float(m_rdf.group(1))
            data['rdf_verdict'] = m_rdf.group(2)

        # --- IMPROVED SYMMETRY PARSING ---
        # 1. Get the Pass/Fail status from the summary
        m_sym_status = re.search(r"Symmetry Retention.*:\s*(PASS|FAIL)", content)
        if m_sym_status:
            data['sym_verdict'] = m_sym_status.group(1)

        # 2. Extract all rows from the table: [Tolerance, Symbol, Number, Match]
        # This regex matches lines like: "1.00    F-43m    216    NO"
        table_matches = re.findall(r"(\d+\.\d+)\s+([A-Za-z0-9/\-_]+)\s+(\d+)\s+(?:YES|NO)", content)
        
        if table_matches:
            # The "Latest" or "Most Stable" entry is the last one (highest tolerance)
            last_entry = table_matches[-1]
            data['sym_group'] = f"{last_entry[1]}"
            
            # Format the full list for the "Sym List" column
            unique_sgs = []
            for _, sg_sym, sg_num in table_matches:
                sg_str = f"{sg_sym} ({sg_num})"
                if sg_str not in unique_sgs:
                    unique_sgs.append(sg_str)
            data["post_md_sgs"] = "; ".join(unique_sgs)
        
    except Exception as e:
        print(f"Error parsing {report_path}: {e}")
        data["status"] = "PARSE_ERROR"
    
    return data

def get_md_info(run_dir):
    md_dir = os.path.join(run_dir, "md_stability_analysis")
    info = {
        'MD Verdict': "N/A", 'MD Confidence': "N/A",
        'RMSD Val': "N/A", 'RMSD Verdict': "N/A",
        'Vol Val': "N/A", 'Vol Verdict': "N/A",
        'RDF Val': "N/A", 'RDF Verdict': "N/A",
        'Sym Verdict': "N/A", 'Sym Group': "N/A", 'Sym List': "N/A",
        'MD Status': "MISSING"
    }
    
    report_glob = glob.glob(os.path.join(md_dir, "*_md_stability_report.txt"))
    if report_glob:
        stats = parse_md_report(report_glob[0])
        info['MD Verdict'] = stats.get("verdict", "N/A")
        info['MD Confidence'] = stats.get("confidence", "N/A")
        
        info['RMSD Val'] = stats.get("rmsd_val", "N/A")
        info['RMSD Verdict'] = stats.get("rmsd_verdict", "N/A")
        
        info['Vol Val'] = stats.get("vol_val", "N/A")
        info['Vol Verdict'] = stats.get("vol_verdict", "N/A")
        
        info['RDF Val'] = stats.get("rdf_val", "N/A")
        info['RDF Verdict'] = stats.get("rdf_verdict", "N/A")
        
        info['Sym Verdict'] = stats.get("sym_verdict", "N/A")
        info['Sym Group'] = stats.get("sym_group", "N/A")
        info['Sym List'] = stats.get("post_md_sgs", "N/A")
        
        info['MD Status'] = stats.get("status", "COMPLETED")
    else:
        if os.path.exists(md_dir):
            info['MD Status'] = "INCOMPLETE"
        else:
            info['MD Status'] = "MISSING"
    
    return info

def check_run_completeness(run_dir):
    issues = []
    relax_dir = os.path.join(run_dir, "initial_relaxation_for_single_run")
    md_dir = os.path.join(run_dir, "md_stability_analysis")
    
    if not os.path.exists(relax_dir):
        issues.append("missing_relax_dir")
    
    md_report_exists = False
    if not os.path.exists(md_dir):
        issues.append("missing_md_dir")
    else:
        report_glob = glob.glob(os.path.join(md_dir, "*_md_stability_report.txt"))
        if report_glob:
            md_report_exists = True

    if not md_report_exists:
        found_logs = []
        found_logs.extend(glob.glob(os.path.join(run_dir, "*.log")))
        engine_dir = os.path.dirname(run_dir.rstrip('/'))
        log_dirs = glob.glob(os.path.join(engine_dir, "batch_logs*"))
        
        folder_name = os.path.basename(run_dir.rstrip('/'))
        run_signature = folder_name
        
        run_signature = re.sub(r'^MD_cifs_vibroml\d*_(?:mace|macemp|esen|uma)_', '', run_signature, flags=re.IGNORECASE)
        run_signature = re.sub(r'^md_cifs_vibroml\d*_(?:mace|macemp|esen|uma)_', '', run_signature, flags=re.IGNORECASE)
        run_signature = re.sub(r'^MD_cifs_', '', run_signature, flags=re.IGNORECASE)
        run_signature = re.sub(r'_(?:mace|macemp|esen|uma)_md_stability.*$', '', run_signature, flags=re.IGNORECASE)
        run_signature = re.sub(r'_(?:mace|macemp|esen|uma)_md.*$', '', run_signature, flags=re.IGNORECASE)
        run_signature = re.sub(r'_\d{8}-\d{6}$', '', run_signature)
        run_signature = re.sub(r'_phonon_output$', '', run_signature)
        
        if run_signature and log_dirs:
            target_sig = run_signature.lower()
            for ld in log_dirs:
                if os.path.exists(ld):
                    try:
                        all_files = os.listdir(ld)
                        for f in all_files:
                            if f.endswith('.log') and (target_sig in f.lower()):
                                found_logs.append(os.path.join(ld, f))
                    except OSError: pass

        if found_logs:
            found_logs.sort(key=os.path.getmtime)
            latest_log = found_logs[-1]
            try:
                with open(latest_log, 'r', errors='replace') as f:
                    log_content = f.read()
                    if "DUE TO TIME LIMIT" in log_content or "CANCELLED AT" in log_content:
                        issues.append("slurm_timeout")
                    elif "No edges found in input system" in log_content:
                        issues.append("graph_break_unstable")
                    elif "torch.OutOfMemoryError" in log_content or "CUDA out of memory" in log_content:
                        issues.append("cuda_oom")
                    elif ("overflow encountered" in log_content or "invalid value encountered" in log_content or "unsupported operand type(s) for *" in log_content):
                        issues.append("md_explosion")
                    else:
                        issues.append("incomplete_no_error_found")
            except Exception:
                issues.append("log_read_error")
        else:
            issues.append("no_logs_found")

    return issues

def process_run(path, engine):
    folder_name = os.path.basename(path.rstrip('/'))
    data = {'Engine': engine, 'Run ID': clean_run_id(folder_name)}
    
    data.update(parse_folder_name_info(folder_name))
    data.update(get_relaxed_info(path))
    data.update(get_md_info(path))
    
    issues = check_run_completeness(path)
    
    rmsd = data.get('RMSD Val')
    if isinstance(rmsd, (int, float)) and rmsd > 100.0:
        if "md_explosion" not in issues: issues.append("md_explosion")
            
    data['Completeness Issues'] = "; ".join(issues) if issues else "COMPLETE"
    
    if data['Completeness Issues'] != "COMPLETE":
        if "graph_break_unstable" in issues or "md_explosion" in issues:
            data['Overall Status'] = "FAILED_UNSTABLE"
        else:
            data['Overall Status'] = "INCOMPLETE"
    elif data.get('MD Status') == "FAILED":
        data['Overall Status'] = "FAILED"
    elif data.get('MD Verdict') == "STABLE":
        data['Overall Status'] = "STABLE"
    elif data.get('MD Verdict') == "UNSTABLE":
        data['Overall Status'] = "UNSTABLE"
    else:
        data['Overall Status'] = "UNKNOWN"
    
    if data['Energy (eV/atom)'] is None:
        data['Energy (eV/atom)'] = data.get('Relaxed Energy (eV/atom)')
    
    if not data['Run ID'] or data['Run ID'] == data.get('Composition', ''):
        comp = data.get('Composition', 'Unknown')
        data['Run ID'] = f"{comp}_MD"
    
    return data

def report_missing_runs(df, expected_engines=None):
    if expected_engines is None:
        expected_engines = ENGINE_SORT_ORDER

    print(f"\n{'='*60}")
    print(f"--- Missing/Incomplete Runs Report ---")
    print(f"Expected Engines: {', '.join(expected_engines)}")
    print(f"{'='*60}")

    incomplete = df[
        df['Overall Status'].isin(['INCOMPLETE', 'FAILED', 'UNKNOWN']) & 
        (~df['Completeness Issues'].str.contains("graph_break_unstable", na=False)) & 
        (~df['Completeness Issues'].str.contains("md_explosion", na=False)) &
        (~df['Completeness Issues'].str.contains("cuda_oom", na=False))
    ]
    
    if incomplete.empty:
        print("✅ All runs are complete or have valid physical failures!")
    else:
        print(f"\n⚠️  Found {len(incomplete)} runs with issues needing attention:")
        print(incomplete[['Engine', 'Run ID', 'Composition', 'Overall Status', 
                          'Completeness Issues', 'MD Verdict']].to_markdown(index=False))
        
        print(f"\n--- Issues by Type ---")
        issue_counts = incomplete['Completeness Issues'].value_counts()
        for issue, count in issue_counts.items():
            print(f"  {issue}: {count}")
    
    graph_break_runs = df[df['Completeness Issues'].str.contains("graph_break_unstable", na=False)]
    if not graph_break_runs.empty:
        print(f"\n--- Graph Break Detected (Physically Unstable) ---")
        print(f"  {len(graph_break_runs)} runs identified as unstable (atoms exceeded model cutoff):")
        for _, row in graph_break_runs.iterrows():
            print(f"    - {row['Engine']}: {row['Run ID']} ({row['Composition']})")

    explosion_runs = df[df['Completeness Issues'].str.contains("md_explosion", na=False)]
    if not explosion_runs.empty:
        print(f"\n💥 --- MD Explosion Detected (Physically Unstable) ---")
        print(f"  {len(explosion_runs)} runs failed due to catastrophic instability (Atoms flying apart):")
        for _, row in explosion_runs.iterrows():
            val = row['RMSD Val'] if row['RMSD Val'] != "N/A" else "Overflow/Crash"
            print(f"    - {row['Engine']}: {row['Run ID']} ({row['Composition']}) - RMSD: {val}")

    oom_runs = df[df['Completeness Issues'].str.contains("cuda_oom", na=False)]
    if not oom_runs.empty:
        print(f"\n--- CUDA OOM Detected (Infrastructure Issue) ---")
        print(f"  {len(oom_runs)} runs failed due to GPU memory exhaustion:")
        for _, row in oom_runs.iterrows():
            print(f"    - {row['Engine']}: {row['Run ID']} ({row['Composition']})")
        print(f"  NOTE: These should be retried with PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True")
    
    csv_out = "vibroml_md_incomplete_runs.csv"
    incomplete.to_csv(csv_out, index=False)
    print(f"\nIncomplete runs report saved to: {csv_out}")
    return incomplete

def report_rerun_candidates(df):
    print(f"\n{'='*60}")
    print(f"--- Rerun Candidates ---")
    print(f"{'='*60}")
    
    needs_rerun = df[
        (df['Completeness Issues'] != "COMPLETE") & 
        (~df['Completeness Issues'].str.contains("graph_break_unstable", na=False)) & 
        (~df['Completeness Issues'].str.contains("md_explosion", na=False)) &
        (~df['Completeness Issues'].str.contains("cuda_oom", na=False))
    ]
    
    if needs_rerun.empty:
        print("✅ No runs need to be rerun!")
        return
    
    print(f"\n🔄 {len(needs_rerun)} runs need to be rerun:")
    for engine in ENGINE_SORT_ORDER:
        engine_runs = needs_rerun[needs_rerun['Engine'] == engine]
        if not engine_runs.empty:
            print(f"\n  {engine}: {len(engine_runs)} runs")
            for _, row in engine_runs.iterrows():
                print(f"    - {row['Run ID']} ({row['Composition']})")
                if row['Completeness Issues'] != "COMPLETE":
                    print(f"      Issues: {row['Completeness Issues']}")
    
    rerun_csv = "vibroml_md_rerun_list.csv"
    needs_rerun[['Engine', 'Run ID', 'Composition', 'Completeness Issues', 'Energy (eV/atom)']].to_csv(rerun_csv, index=False)
    print(f"\nRerun list saved to: {rerun_csv}")
    
    explosion_runs = df[df['Completeness Issues'].str.contains("md_explosion", na=False)]
    if not explosion_runs.empty:
        print(f"\n💥 {len(explosion_runs)} runs skipped (md_explosion - Physical Instability).")

    timeout_runs = df[df['Completeness Issues'].str.contains("slurm_timeout", na=False)]
    if not timeout_runs.empty:
        print(f"\n⏳ {len(timeout_runs)} runs timed out (Slurm Time Limit):")
        for _, row in timeout_runs.iterrows():
            print(f"    - {row['Engine']}: {row['Run ID']} ({row['Composition']})")
        print(f"  ADVICE: Increase the walltime in your SBATCH script for these runs.")
        
    return needs_rerun

def main():
    all_data = []
    cwd = os.getcwd()
    
    expected_by_parent = {}
    if os.path.exists(ROOT_CIF_DIR):
        for root, dirs, files in os.walk(ROOT_CIF_DIR):
            parent = os.path.basename(root)
            for f in files:
                if f.endswith('.cif'):
                    if parent not in expected_by_parent:
                        expected_by_parent[parent] = []
                    expected_by_parent[parent].append(f.replace('.cif', ''))
    
    total_cifs = sum(len(v) for v in expected_by_parent.values())
    print(f"Scanning directories in {cwd}...")
    print(f"Looking for MD evaluation directories: {list(ENGINE_DIRS.keys())}")
    if expected_by_parent:
        print(f"Reference CIFs found: {total_cifs} in {ROOT_CIF_DIR}")
    
    for dir_name, engine_name in ENGINE_DIRS.items():
        full_dir_path = os.path.join(cwd, dir_name)
        if os.path.exists(full_dir_path):
            print(f"Processing {engine_name} MD runs in {dir_name}...")
            processed = 0
            for run_folder in os.listdir(full_dir_path):
                run_path = os.path.join(full_dir_path, run_folder)
                if os.path.isdir(run_path) and run_folder.startswith("MD_cifs"):
                    try:
                        data = process_run(run_path, engine_name)
                        all_data.append(data)
                        processed += 1
                    except Exception as e:
                        print(f"  Error reading {run_folder}: {e}")
            print(f"  Found {processed} potential run folders in {engine_name}")
        else:
            print(f"Skipping {dir_name} (not found)")

    if not all_data:
        print("No MD data found. Exiting.")
        return

    # --- DEDUPLICATION LOGIC ---
    df = pd.DataFrame(all_data)
    
    def rank_quality(row):
        issues = row['Completeness Issues']
        if issues == "COMPLETE": return 0
        elif "graph_break_unstable" in issues or "md_explosion" in issues: return 1
        elif "slurm_timeout" in issues or "cuda_oom" in issues: return 2
        else: return 3

    print(f"\n{'='*60}")
    print(f"--- Handling Duplicates ---")
    print(f"{'='*60}")
    
    initial_count = len(df)
    df['Quality_Rank'] = df.apply(rank_quality, axis=1)
    df = df.sort_values(by=['Engine', 'Run ID', 'Quality_Rank', 'Energy (eV/atom)'], ascending=[True, True, True, True])
    df_dedup = df.drop_duplicates(subset=['Engine', 'Run ID'], keep='first')
    
    print(f"Entries removed: {initial_count - len(df_dedup)}")
    print(f"Unique structures: {len(df_dedup)}")
    df = df_dedup

    # --- NORMALIZE COMPOSITION NAMES TO MATCH TARGET LIST ---
    # This fixes the issue where "FLi" != "LiF" and gets dropped by pd.Categorical
    print("Normalizing composition names...")
    
    # Pre-compute target element sets for faster lookup
    target_comps = {}
    for target in COMPOUND_ORDER:
        try:
            f = Formula(target)
            target_comps[target] = f.count()
        except: pass
        
    def normalize_name(found_name):
        if found_name in COMPOUND_ORDER:
            return found_name
        
        # Check against targets
        try:
            f_found = Formula(found_name)
            c_found = f_found.count()
            
            for target_name, c_target in target_comps.items():
                if c_found == c_target:
                    return target_name
        except: pass
        
        return found_name

    df['Composition'] = df['Composition'].apply(normalize_name)

    # --- SORTING LOGIC FOR FINAL DISPLAY ---
    df['Engine'] = pd.Categorical(df['Engine'], categories=ENGINE_SORT_ORDER, ordered=True)
    df['Composition'] = pd.Categorical(df['Composition'], categories=COMPOUND_ORDER, ordered=True)

    structure_energies = df.groupby('Run ID')['Energy (eV/atom)'].mean().reset_index()
    structure_energies.rename(columns={'Energy (eV/atom)': 'Structure_Mean_Energy'}, inplace=True)
    df = df.merge(structure_energies, on='Run ID', how='left')

    df = df.sort_values(by=['Composition', 'Structure_Mean_Energy', 'Engine'], ascending=[True, False, True])
    
    final_cols = [
        'Engine', 'Run ID', 'Composition', 'Initial Space Group',
        'Energy (eV/atom)', 'Structure_Mean_Energy',
        'Overall Status', 'MD Verdict', 
        'RMSD Val', 'RMSD Verdict', 
        'Vol Val', 'Vol Verdict',
        'RDF Val', 'RDF Verdict',
        'Sym Verdict', 'Sym Group', 'Sym List',
        'Completeness Issues'
    ]
    df_final = df[[c for c in final_cols if c in df.columns]]

    print(f"\n{'='*130}")
    print(f"--- Structure-Centric MD Summary ---")
    print(f"Order: {', '.join(COMPOUND_ORDER)}")
    print(f"Sorting: Higher Energy -> Lower Energy")
    print(f"Format: Value (Verdict) | P=PASS, F=FAIL")
    print(f"{'='*130}")
    
    for comp in COMPOUND_ORDER:
        comp_df = df_final[df_final['Composition'] == comp]
        if comp_df.empty: continue
            
        print(f"\n### {comp} ###")
        unique_runs = comp_df['Run ID'].unique()
        
        for run_id in unique_runs:
            run_data = comp_df[comp_df['Run ID'] == run_id]
            avg_energy = run_data['Structure_Mean_Energy'].iloc[0]
            header_energy = f"{avg_energy:.3f}" if pd.notnull(avg_energy) else "N/A"
            
            print(f"\n  Structure: {run_id:<40} (Avg E: {header_energy} eV/atom)")
            print(f"  {'-'*130}")
            print(f"  {'Engine':<12} | {'Energy':<8} | {'RMSD (A)':<14} | {'Vol (%)':<14} | {'RDF':<14} | {'Sym':<10} | {'Status'}")
            print(f"  {'-'*130}")
            
            for engine in ENGINE_SORT_ORDER:
                row = run_data[run_data['Engine'] == engine]
                if not row.empty:
                    r = row.iloc[0]
                    e_val = f"{r['Energy (eV/atom)']:.3f}" if pd.notnull(r['Energy (eV/atom)']) else "N/A"
                    
                    def fmt_metric(val, verdict):
                        if verdict == "PASS": v_code = "P"
                        elif verdict == "FAIL": v_code = "F"
                        else: v_code = "?"
                        if isinstance(val, (float, int)): return f"{val:.2f} ({v_code})"
                        return f"{str(val)[:5]} ({v_code})"

                    rmsd_str = fmt_metric(r['RMSD Val'], r['RMSD Verdict'])
                    vol_str = fmt_metric(r['Vol Val'], r['Vol Verdict'])
                    rdf_str = fmt_metric(r['RDF Val'], r['RDF Verdict'])
                    sym_str = r['Sym Verdict'] if r['Sym Verdict'] != "N/A" else "-"
                    status = r['Overall Status']
                    if status == "FAILED_UNSTABLE": status = "PHYS_FAIL"
                    
                    print(f"  {engine:<12} | {e_val:<8} | {rmsd_str:<14} | {vol_str:<14} | {rdf_str:<14} | {sym_str:<10} | {status}")
                else:
                    print(f"  {engine:<12} | {'-':<8} | {'-':<14} | {'-':<14} | {'-':<14} | {'-':<10} | NOT_FOUND")
            print(f"  {'-'*130}")
    # Delete the Structure Mean Energy as it is redundant
    df_final = df_final.drop(columns=['Structure_Mean_Energy'])
    csv_name = "vibroml_master_md_summary.csv"
    df_final.to_csv(csv_name, index=False)
    print(f"\nFull results saved to: {csv_name}")

    report_missing_runs(df_final)
    report_rerun_candidates(df_final)

if __name__ == "__main__":
    main()