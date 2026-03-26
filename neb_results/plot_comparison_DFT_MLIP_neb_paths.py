#!/usr/bin/env python3
"""
Generate comparison plots for Phase Transition NEB paths.
ALIGNMENT: Image 0 set to 0.0 eV.
METRICS: 
  - Energy in eV/atom.
  - Forces in eV/A.
  - NEW: "Endpts F Err" = |F_init_DFT - F_init_MLIP| + |F_final_DFT - F_final_MLIP|
  - Ranking tables included in summary.

OUTPUTS:
  1. PNG Plots.
  2. 'all_metrics_summary.txt' with rankings.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
BASE_DIR = "."
VASP_SUBDIR = "vasp_plots_and_neb_data"
OUTPUT_DIR = "comparison_plots"
SUMMARY_FILENAME = "all_metrics_summary.txt"

ATOM_COUNTS = {
    'Bi2Sn2O7': 44,
    'CsPbI3': 20,
    'HfO2': 12,
    'LiF': 4
}

MLIPS = ['MACE', 'ESEN', 'UMA']
COLORS = {'MACE': '#d62728', 'ESEN': '#2ca02c', 'UMA': '#9467bd', 'DFT': '#1f77b4'}
MARKERS = {'MACE': 's', 'ESEN': '^', 'UMA': 'D', 'DFT': 'o'}

# --- Font Sizes ---
FS_TITLE = 24
FS_LABEL = 20
FS_TICK = 18
FS_LEGEND = 18
FS_TEXT = 15

def parse_neb_summary(filepath):
    data = [] 
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        section_found = False
        for line in lines:
            if 'Energy and Force Profile Along Path' in line:
                section_found = True; continue
            if section_found:
                if line.strip().startswith('-') or 'Image' in line or not line.strip(): continue
                parts = line.split()
                if len(parts) >= 5 and parts[0].isdigit():
                    try:
                        data.append({
                            'image': int(parts[0]),
                            'rel_energy_total_ev': float(parts[2]),
                            'max_force': float(parts[4])
                        })
                    except ValueError: continue
    except: return None
    
    data.sort(key=lambda x: x['image'])
    if not data: return None
    return {
        'image': [d['image'] for d in data],
        'rel_energy_total_ev': [d['rel_energy_total_ev'] for d in data],
        'max_force': [d['max_force'] for d in data]
    }

def parse_vasp_csv(filepath):
    data = []
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                data.append({
                    'image': int(row['Image_Index']),
                    'rel_energy_per_atom': float(row['Relative_Energy_eV_per_atom']),
                    'max_force': float(row['Max_Force_eV_A'])
                })
    except: return None

    data.sort(key=lambda x: x['image'])
    if not data: return None
    return {
        'image': [d['image'] for d in data],
        'rel_energy_per_atom': [d['rel_energy_per_atom'] for d in data],
        'max_force': [d['max_force'] for d in data]
    }

def get_aligned_data(structure, mlip_data, vasp_data):
    n_atoms = ATOM_COUNTS.get(structure, 1)
    
    # Process VASP (Already eV/atom)
    v_e = np.array(vasp_data['rel_energy_per_atom'])
    v_e = v_e - v_e[0] 
    
    # Process MLIP (Total eV -> eV/atom)
    m_e = np.array(mlip_data['rel_energy_total_ev']) / n_atoms
    m_e = m_e - m_e[0]
    
    # Process Forces (Intensive, no scaling)
    v_f = np.array(vasp_data['max_force'])
    m_f = np.array(mlip_data['max_force'])
    
    return {
        'images': np.array(vasp_data['image']),
        'vasp_e': v_e, 'vasp_f': v_f,
        'mlip_e': m_e, 'mlip_f': m_f
    }

def calculate_limits(data_cache):
    limits = {}
    for struct in ATOM_COUNTS.keys():
        e_vals, f_vals = [], []
        if struct in data_cache:
            for mlip in data_cache[struct]:
                d = data_cache[struct][mlip]
                e_vals.extend(d['vasp_e']); e_vals.extend(d['mlip_e'])
                f_vals.extend(d['vasp_f']); f_vals.extend(d['mlip_f'])
        if e_vals:
            e_min, e_max = min(e_vals), max(e_vals)
            f_min, f_max = min(f_vals), max(f_vals)
            e_rng = e_max - e_min if e_max != e_min else 0.1
            f_rng = f_max - f_min if f_max != f_min else 0.1
            
            limits[struct] = {
                'e_lim': (e_min - (e_rng*0.05), e_max + (e_rng*0.35)),
                'f_lim': (f_min - (f_rng*0.05), f_max + (f_rng*0.35))
            }
    return limits

def generate_report_entry(structure, mlip, images, me, ve, mf, vf, metrics):
    lines = []
    lines.append("-" * 80)
    lines.append(f"Structure: {structure}  |  Model: {mlip}")
    lines.append("-" * 80)
    lines.append("SCALAR METRICS:")
    lines.append(f"  Reaction Energy Error:  {metrics['err_rxn']:.6f} eV/atom")
    lines.append(f"  Barrier Energy Error:   {metrics['err_barrier']:.6f} eV/atom")
    lines.append(f"  DFT Barrier Height:     {metrics['barrier_dft']:.6f} eV/atom")
    lines.append(f"  MAE Energy:             {metrics['mae_e']:.6f} eV/atom")
    lines.append(f"  Avg Path Energy (DFT):  {metrics['avg_path_dft']:.6f} eV/atom")
    lines.append(f"  MAE Force:              {metrics['mae_f']:.6f} eV/Å")
    lines.append(f"  Endpts Force Err:       {metrics['err_endpts_f']:.6f} eV/Å (Init+Final)")
    lines.append("")
    
    lines.append("DETAILED DATA (Aligned):")
    headers = ["Img", "E_MLIP", "E_DFT", "E_Err", "F_MLIP", "F_DFT", "F_Err"]
    lines.append(f"{headers[0]:<4} {headers[1]:<10} {headers[2]:<10} {headers[3]:<10} "
                 f"{headers[4]:<10} {headers[5]:<10} {headers[6]:<10}")
    
    for i in range(len(images)):
        e_err = abs(me[i] - ve[i])
        f_err = abs(mf[i] - vf[i])
        lines.append(
            f"{images[i]:<4d} {me[i]:<10.5f} {ve[i]:<10.5f} {e_err:<10.5f} "
            f"{mf[i]:<10.5f} {vf[i]:<10.5f} {f_err:<10.5f}"
        )
    lines.append("\n")
    return "\n".join(lines)

def generate_ranking_block(structure, structure_metrics):
    """Generates a sorted ranking table for the structure."""
    if not structure_metrics: return ""
    
    lines = []
    lines.append("=" * 80)
    lines.append(f"PERFORMANCE RANKING FOR {structure}")
    lines.append("=" * 80)
    
    # 1. Path Quality (DFT Path Avg)
    sorted_path = sorted(structure_metrics.items(), key=lambda x: x[1]['avg_path_dft'])
    lines.append("By Path Quality (DFT Avg Energy - Lower is Better):")
    for i, (m, v) in enumerate(sorted_path):
        lines.append(f"  {i+1}. {m:<6} : {v['avg_path_dft']:.5f} eV/at (Barrier: {v['barrier_dft']:.4f})")
    
    # 2. Accuracy (MAE Energy)
    sorted_mae_e = sorted(structure_metrics.items(), key=lambda x: x[1]['mae_e'])
    lines.append("\nBy Accuracy (MAE Energy):")
    for i, (m, v) in enumerate(sorted_mae_e):
        lines.append(f"  {i+1}. {m:<6} : {v['mae_e']:.5f} eV/at")

    # 3. Stability (Endpoints Force Error)
    sorted_stab = sorted(structure_metrics.items(), key=lambda x: x[1]['err_endpts_f'])
    lines.append("\nBy Stability (Endpoint Forces Error):")
    for i, (m, v) in enumerate(sorted_stab):
        lines.append(f"  {i+1}. {m:<6} : {v['err_endpts_f']:.5f} eV/Å")

    # 4. Accuracy (MAE Force)
    sorted_mae_f = sorted(structure_metrics.items(), key=lambda x: x[1]['mae_f'])
    lines.append("\nBy Accuracy (MAE Force):")
    for i, (m, v) in enumerate(sorted_mae_f):
        lines.append(f"  {i+1}. {m:<6} : {v['mae_f']:.5f} eV/Å")

    lines.append("\n" + "="*80 + "\n\n")
    return "\n".join(lines)

def plot_comparison(structure, mlip, data, limits, output_dir):
    images = data['images']
    ve, vf = data['vasp_e'], data['vasp_f']
    me, mf = data['mlip_e'], data['mlip_f']
    
    # --- Metrics ---
    err_rxn = abs(ve[-1] - me[-1])
    err_barrier = abs(max(ve) - max(me))
    
    # Endpoint Stability: Sum of errors at Image 0 and Image N
    err_init_f = abs(vf[0] - mf[0])
    err_term_f = abs(vf[-1] - mf[-1])
    err_endpts_f = err_init_f + err_term_f
    
    n = min(len(ve), len(me))
    mae_e = np.mean(np.abs(ve[1:n] - me[1:n])) if n > 1 else 0.0
    mae_f = np.mean(np.abs(vf[:n] - mf[:n]))

    avg_path_dft = np.mean(ve)
    avg_path_mlip = np.mean(me)
    barrier_dft = max(ve) 

    metrics_dict = {
        'err_rxn': err_rxn,
        'err_barrier': err_barrier,
        'err_endpts_f': err_endpts_f,
        'mae_e': mae_e,
        'mae_f': mae_f,
        'avg_path_dft': avg_path_dft,
        'barrier_dft': barrier_dft
    }
    
    report_text = generate_report_entry(structure, mlip, images, me, ve, mf, vf, metrics_dict)

    # --- Plotting ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Energy
    ax1.plot(images, ve, marker=MARKERS['DFT'], color=COLORS['DFT'], lw=4, label='DFT', ms=12)
    ax1.plot(images, me, marker=MARKERS.get(mlip, 's'), color=COLORS.get(mlip, 'r'), ls='--', lw=4, label=mlip, ms=12)
    
    ax1.set_xlabel("Image Index", fontsize=FS_LABEL)
    ax1.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
    ax1.set_title(f"{structure} - Energy", fontsize=FS_TITLE)
    ax1.tick_params(axis='both', which='major', labelsize=FS_TICK)
    ax1.legend(fontsize=FS_LEGEND, loc='upper left')
    ax1.grid(True, alpha=0.3); ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
    if structure in limits: ax1.set_ylim(limits[structure]['e_lim'])

    # Energy Text
    lines_e = [
        f"Rxn E Err:    {err_rxn:.4f}",
        f"Barrier Err:  {err_barrier:.4f}",
        f"MAE Energy:   {mae_e:.4f}",
        f"DFT Path Avg: {avg_path_dft:.4f}"
    ]
    if barrier_dft > 0.001:
        lines_e.append(f"DFT Barrier:  {barrier_dft:.4f}")
        
    stats_e = "\n".join(lines_e)
    
    ax1.text(0.96, 0.96, stats_e, transform=ax1.transAxes, fontsize=FS_TEXT,
             verticalalignment='top', horizontalalignment='right', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='#aaaaaa'))

    # Force
    ax2.plot(images, vf, marker=MARKERS['DFT'], color=COLORS['DFT'], lw=4, label='DFT', ms=12)
    ax2.plot(images, mf, marker=MARKERS.get(mlip, 's'), color=COLORS.get(mlip, 'r'), ls='--', lw=4, label=mlip, ms=12)
    
    ax2.set_xlabel("Image Index", fontsize=FS_LABEL)
    ax2.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
    ax2.set_title(f"{structure} - Forces", fontsize=FS_TITLE)
    ax2.tick_params(axis='both', which='major', labelsize=FS_TICK)
    ax2.legend(fontsize=FS_LEGEND, loc='upper left') 
    ax2.grid(True, alpha=0.3)
    if structure in limits: ax2.set_ylim(limits[structure]['f_lim'])

    stats_f = (f"Endpts F Err: {err_endpts_f:.4f}\n"
               f"MAE Force:    {mae_f:.4f}")
               
    ax2.text(0.96, 0.96, stats_f, transform=ax2.transAxes, fontsize=FS_TEXT,
             verticalalignment='top', horizontalalignment='right', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='#aaaaaa'))

    plt.tight_layout()
    out_name = f"{structure}_{mlip}_comparison.png"
    plt.savefig(os.path.join(output_dir, out_name), dpi=150)
    plt.close()
    print(f"Created: {out_name}")
    
    return report_text, metrics_dict

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Loading data...")
    data_cache = {} 
    
    for struct in ATOM_COUNTS:
        data_cache[struct] = {}
        for mlip in MLIPS:
            neb_file = f"{struct}_{mlip}_neb_summary.txt"
            vasp_file = os.path.join(VASP_SUBDIR, f"{struct}_{mlip}_NEB_plot_data.csv")
            m, v = parse_neb_summary(neb_file), parse_vasp_csv(vasp_file)
            if m and v and len(m['image']) == len(v['image']):
                data_cache[struct][mlip] = get_aligned_data(struct, m, v)

    limits = calculate_limits(data_cache)
    full_report_content = []
    
    print("-" * 60)
    
    for struct in sorted(data_cache.keys()):
        struct_metrics = {}
        
        # 1. Process Individual MLIPs
        for mlip in MLIPS: 
            if mlip in data_cache[struct]:
                txt, metrics = plot_comparison(struct, mlip, data_cache[struct][mlip], limits, OUTPUT_DIR)
                full_report_content.append(txt)
                struct_metrics[mlip] = metrics
        
        # 2. Generate Ranking for this Structure
        ranking_txt = generate_ranking_block(struct, struct_metrics)
        full_report_content.append(ranking_txt)

    summary_path = os.path.join(OUTPUT_DIR, SUMMARY_FILENAME)
    with open(summary_path, 'w') as f:
        f.write("\n".join(full_report_content))
        
    print("-" * 60)
    print(f"Summary saved to: {summary_path}")

if __name__ == "__main__":
    main()