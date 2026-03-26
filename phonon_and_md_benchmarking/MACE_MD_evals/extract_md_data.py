import os
import re
import glob
import pandas as pd
from ase.io import read
from math import gcd
from functools import reduce

def get_minimal_formula(cif_path):
    """
    Reads a CIF and returns the minimal (empirical) Hill formula.
    Example: Li64F64 -> FLi (Alphabetical sort for consistency)
    """
    if not os.path.exists(cif_path):
        return "N/A"
    try:
        atoms = read(cif_path)
        symbols = atoms.get_chemical_symbols()
        counts = {}
        for s in symbols:
            counts[s] = counts.get(s, 0) + 1
        
        if not counts: return "Empty"
        
        common_divisor = reduce(gcd, counts.values())
        reduced_counts = {k: v // common_divisor for k, v in counts.items()}
        
        # Sort alphabetically
        sorted_elements = sorted(reduced_counts.keys())
        
        formula = ""
        for el in sorted_elements:
            n = reduced_counts[el]
            formula += f"{el}"
            if n > 1: formula += f"{n}"
        return formula
    except Exception:
        return "Error"

def parse_detailed_symmetry_file(filepath):
    """Parses multi-line symmetry analysis text file."""
    if not os.path.exists(filepath): return "N/A"
    try:
        with open(filepath, 'r') as f: content = f.read()
        
        sym_match = re.search(r"International symbol:\s+(.+)", content)
        num_match = re.search(r"Space group number:\s+(\d+)", content)
        
        if sym_match and num_match:
            return f"{sym_match.group(1).strip()} ({num_match.group(1)})"
        elif sym_match:
            return sym_match.group(1).strip()
        return "Unknown"
    except Exception: return "Error"

def clean_run_id(folder_name):
    """
    Extracts a concise Run ID from the folder name.
    Prioritizes: unique/top identifiers > sample/iter > stripped name.
    """
    # 1. Look for specific pattern: unique_#_iter#_sample# OR top_#_iter#...
    # Matches: unique_6_iter1_sample16, top_1_iter2_sample4, targeted_GA_top_1
    specific_id_match = re.search(r'((?:OPTRANDOM_)?(?:unique|top|targeted_GA)_\d+(?:_sg_[A-Za-z0-9]+)?(?:_iter\d+)?(?:_sample\d+)?)', folder_name)
    
    if specific_id_match:
        return specific_id_match.group(1)

    # 2. Fallback: Heuristic cleaning
    # Remove standard prefixes
    clean = re.sub(r'^MD_cifs_vibroml\d*_mace_', '', folder_name)
    clean = re.sub(r'^MD_cifs_', '', clean)
    
    # Remove standard suffixes (timestamps, phonon output labels)
    clean = re.sub(r'_MACE_MD_STABILITY.*$', '', clean)
    clean = re.sub(r'_MACE_MD.*$', '', clean) # shorter variant
    clean = re.sub(r'_\d{8}-\d{6}$', '', clean) # timestamp only
    
    return clean

def extract_energy_from_folder(folder_name):
    """Extracts energy from folder string even if 'energy_' prefix is missing."""
    # Matches _m6p031730 or _p4p8197 preceded by underscore or 'energy_'
    match = re.search(r'(?:_|energy_)([mp])(\d+)p(\d+)', folder_name)
    if match:
        sign_char = match.group(1)
        integer_part = match.group(2)
        decimal_part = match.group(3)
        value = float(f"{integer_part}.{decimal_part}")
        if sign_char == 'm': value = -value
        return value
    return None

def parse_md_report(report_path):
    """Parses MD report for verdicts and metrics."""
    data = {
        "verdict": "N/A", "confidence": "N/A",
        "rmsd_verdict": "N/A", "rmsd_val": "N/A",
        "vol_verdict": "N/A", "vol_val": "N/A",
        "rdf_verdict": "N/A", "rdf_val": "N/A",
        "sym_verdict": "N/A", "post_md_sgs": []
    }
    
    if not os.path.exists(report_path): return data

    with open(report_path, 'r') as f: content = f.read()

    # Verdicts
    v_match = re.search(r"Stability Verdict:\s+(\w+)", content)
    c_match = re.search(r"Confidence Level:\s+(\w+)", content)
    data["verdict"] = v_match.group(1) if v_match else "N/A"
    data["confidence"] = c_match.group(1) if c_match else "N/A"

    # Criteria & Metrics
    # RMSD
    rmsd_val = re.search(r"Mean RMSD:\s+([\d\.]+)\s+Å", content)
    if rmsd_val: data["rmsd_val"] = float(rmsd_val.group(1))
    
    if re.search(r"RMSD:.*-\s*PASS", content): data["rmsd_verdict"] = "PASS"
    elif re.search(r"RMSD:.*-\s*FAIL", content): data["rmsd_verdict"] = "FAIL"

    # Volume
    vol_val = re.search(r"Max volume fluctuation:\s+([\d\.]+)\%", content)
    if vol_val: data["vol_val"] = float(vol_val.group(1))
    
    if re.search(r"Volume change:.*-\s*PASS", content): data["vol_verdict"] = "PASS"
    elif re.search(r"Volume change:.*-\s*FAIL", content): data["vol_verdict"] = "FAIL"

    # RDF
    rdf_val = re.search(r"Initial vs final correlation:\s+([\d\.]+)", content)
    if rdf_val: data["rdf_val"] = float(rdf_val.group(1))

    if re.search(r"RDF correlation:.*-\s*PASS", content): data["rdf_verdict"] = "PASS"
    elif re.search(r"RDF correlation:.*-\s*FAIL", content): data["rdf_verdict"] = "FAIL"

    # Symmetry
    sym_match = re.search(r"Symmetry Retention.*:\s+(\w+)", content)
    if sym_match:
        data["sym_verdict"] = "PASS" if sym_match.group(1) == "RETAINED" else "FAIL"

    # Found SGs
    found_sgs = re.findall(r"\d+\.\d+\s+([A-Za-z0-9/\-]+)\s+(\d+)\s+(?:YES|NO)", content)
    unique_sgs = []
    for sg_sym, sg_num in found_sgs:
        sg_str = f"{sg_sym} ({sg_num})"
        if sg_str not in unique_sgs and sg_sym != "FAILED":
            unique_sgs.append(sg_str)
    data["post_md_sgs"] = unique_sgs
    return data

def main():
    base_path = "."
    run_folders = [d for d in os.listdir(base_path) if os.path.isdir(d) and d.startswith("MD_cifs")]
    
    results = []
    print(f"Found {len(run_folders)} run directories. Processing...")

    for folder in sorted(run_folders):
        row = {}
        
        # 1. Clean Run ID & Energy
        row["Run_ID"] = clean_run_id(folder)
        row["Energy (eV/atom)"] = extract_energy_from_folder(folder)

        # Paths
        relax_dir = os.path.join(folder, "initial_relaxation_for_single_run")
        md_dir = os.path.join(folder, "md_stability_analysis")

        # 2. Composition (Minimal)
        cif_files = glob.glob(os.path.join(relax_dir, "*_relaxed.cif"))
        row["Composition"] = get_minimal_formula(cif_files[0]) if cif_files else "N/A"

        # 3. Symmetry
        row["Initial_SG"] = parse_detailed_symmetry_file(os.path.join(relax_dir, "initial_symmetry_analysis.txt"))
        row["Relaxed_SG"] = parse_detailed_symmetry_file(os.path.join(relax_dir, "relaxed_symmetry_analysis.txt"))

        # 4. MD Stability
        report_glob = glob.glob(os.path.join(md_dir, "*_md_stability_report.txt"))
        if report_glob:
            stats = parse_md_report(report_glob[0])
            row["Verdict"] = stats["verdict"]
            row["Confidence"] = stats["confidence"]
            
            row["RMSD_Status"] = stats["rmsd_verdict"]
            row["RMSD_Mean(Å)"] = stats["rmsd_val"]
            
            row["Vol_Status"] = stats["vol_verdict"]
            row["Vol_Max_Change(%)"] = stats["vol_val"]
            
            row["RDF_Status"] = stats["rdf_verdict"]
            row["RDF_Corr"] = stats["rdf_val"]
            
            row["Sym_Status"] = stats["sym_verdict"]
            row["MD_Found_SGs"] = "; ".join(stats["post_md_sgs"])
        else:
            for k in ["Verdict", "Confidence", "RMSD_Status", "RMSD_Mean(Å)", "Vol_Status", 
                      "Vol_Max_Change(%)", "RDF_Status", "RDF_Corr", "Sym_Status", "MD_Found_SGs"]:
                row[k] = "N/A"

        results.append(row)

    # DataFrame creation
    df = pd.DataFrame(results)
    
    # Sort: Composition (A-Z) -> Energy (Low-High)
    df.sort_values(by=["Composition", "Energy (eV/atom)"], ascending=[True, True], inplace=True)

    # Reorder columns
    cols = [
        "Run_ID", "Composition", "Energy (eV/atom)", 
        "Initial_SG", "Relaxed_SG", 
        "Verdict", "Confidence", 
        "RMSD_Status", "RMSD_Mean(Å)", 
        "Vol_Status", "Vol_Max_Change(%)", 
        "RDF_Status", "RDF_Corr", 
        "Sym_Status", "MD_Found_SGs"
    ]
    df = df[[c for c in cols if c in df.columns]]

    csv_name = "vibroml_md_summary.csv"
    df.to_csv(csv_name, index=False)
    
    print(f"\nSuccessfully processed {len(df)} runs.")
    print(f"Summary saved to: {csv_name}")
    
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 1000)
    print("\nPreview:")
    print(df.head(10))

if __name__ == "__main__":
    main()