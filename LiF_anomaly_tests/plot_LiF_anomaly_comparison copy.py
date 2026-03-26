#!/usr/bin/env python3
"""
Generate comparison plots for LiF Phase Transition NEB paths (P63mc to Cm, P1, Pm).
ALIGNMENT: Image 0 set to 0.0 eV.
METRICS: 
  - Energy in eV/atom.
  - Forces in eV/A.
  - "Endpts F Err" = |F_init_DFT - F_init_MLIP| + |F_final_DFT - F_final_MLIP|
  - Ranking tables included in summary.

OUTPUTS:
  1. PNG and SVG Plots for each path and MLIP comparison.
  2. 'all_metrics_summary.txt' with rankings.
  3. Combined overview plot with all MLIPs for each path.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt

# --- Configuration ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "comparison_plots")
SUMMARY_FILENAME = "all_metrics_summary.txt"

# Paths (polymorphs) to analyze
PATHS = ['LiF_Cm', 'LiF_P1', 'LiF_Pm']

# Mapping path folder names to LaTeX phase labels for plots
PHASE_LABELS = {
    'LiF_Cm': r'$Cm$',
    'LiF_P1': r'$P1$',
    'LiF_Pm': r'$Pm$'
}
START_PHASE_LABEL = r'$P6_3mc$'

# MLIP models to compare
MLIPS = ['MACE-OMAT', 'MACE-MP', 'ESEN', 'UMA']

# Mapping from MLIP name to result file
MLIP_FILES = {
    'MACE-OMAT': 'results_mace_omat.csv',
    'MACE-MP': 'results_mace_mp.csv',
    'ESEN': 'results_esen.csv',
    'UMA': 'results_uma.csv'
}

# Colors and markers for each MLIP
COLORS = {
    'MACE-OMAT': '#d62728',  # Red
    'MACE-MP': '#ff7f0e',    # Orange
    'ESEN': '#2ca02c',       # Green
    'UMA': '#9467bd',        # Purple
    'DFT': '#1f77b4'         # Blue
}

MARKERS = {
    'MACE-OMAT': 's',
    'MACE-MP': 'p',
    'ESEN': '^',
    'UMA': 'D',
    'DFT': 'o'
}

# --- Font Sizes ---
FS_TITLE = 24
FS_LABEL = 20
FS_TICK = 18
FS_LEGEND = 18
FS_TEXT = 15


def parse_dft_csv(filepath):
    """Parse the DFT results CSV file."""
    data = {}
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                group = row['Group']
                image_name = row['Image']
                # Extract image number from filename
                image_num = int(image_name.split('_')[-1])
                
                if group not in data:
                    data[group] = []
                
                data[group].append({
                    'image': image_num,
                    'total_energy': float(row['Total_Energy_eV']),
                    'energy_per_atom': float(row['Energy_per_atom_eV']),
                    'max_force': float(row['Max_Force_eV_A']),
                    'num_atoms': int(row['Num_Atoms'])
                })
        
        # Sort each group by image number
        for group in data:
            data[group].sort(key=lambda x: x['image'])
        
        return data
    except Exception as e:
        print(f"Error parsing DFT CSV: {e}")
        return None


def parse_mlip_csv(filepath, energy_col, force_col):
    """Parse MLIP results CSV file."""
    data = []
    if not os.path.exists(filepath):
        return None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                filename = row['filename']
                # Extract image number from filename
                image_num = int(filename.split('_')[-1].replace('.cif', ''))
                
                data.append({
                    'image': image_num,
                    'energy_per_atom': float(row[energy_col]),
                    'max_force': float(row[force_col])
                })
        
        # Sort by image number
        data.sort(key=lambda x: x['image'])
        return data
    except Exception as e:
        print(f"Error parsing MLIP CSV {filepath}: {e}")
        return None


def get_aligned_data(dft_data, mlip_data):
    """
    Align DFT and MLIP data.
    Returns aligned energies (relative to image 0) and forces.
    """
    # Extract arrays
    dft_e = np.array([d['energy_per_atom'] for d in dft_data])
    dft_f = np.array([d['max_force'] for d in dft_data])
    images = np.array([d['image'] for d in dft_data])
    
    mlip_e = np.array([d['energy_per_atom'] for d in mlip_data])
    mlip_f = np.array([d['max_force'] for d in mlip_data])
    
    # Align energies to image 0
    dft_e_aligned = dft_e - dft_e[0]
    mlip_e_aligned = mlip_e - mlip_e[0]
    
    return {
        'images': images,
        'dft_e': dft_e_aligned,
        'dft_f': dft_f,
        'mlip_e': mlip_e_aligned,
        'mlip_f': mlip_f
    }


def calculate_limits(data_cache):
    """Calculate y-axis limits for consistent plotting."""
    limits = {}
    for path in PATHS:
        e_vals, f_vals = [], []
        if path in data_cache:
            for mlip in data_cache[path]:
                d = data_cache[path][mlip]
                e_vals.extend(d['dft_e'])
                e_vals.extend(d['mlip_e'])
                f_vals.extend(d['dft_f'])
                f_vals.extend(d['mlip_f'])
        
        if e_vals:
            e_min, e_max = min(e_vals), max(e_vals)
            f_min, f_max = min(f_vals), max(f_vals)
            e_rng = e_max - e_min if e_max != e_min else 0.1
            f_rng = f_max - f_min if f_max != f_min else 0.1
            
            limits[path] = {
                'e_lim': (e_min - (e_rng * 0.05), e_max + (e_rng * 0.35)),
                'f_lim': (f_min - (f_rng * 0.05), f_max + (f_rng * 0.35))
            }
    return limits


def generate_report_entry(path, mlip, images, me, ve, mf, vf, metrics):
    """Generate a text report entry for a path/MLIP combination."""
    lines = []
    lines.append("-" * 80)
    lines.append(f"Path: {path} ({START_PHASE_LABEL} -> {PHASE_LABELS.get(path, path)})  |  Model: {mlip}")
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


def generate_ranking_block(path, path_metrics):
    """Generates a sorted ranking table for the path."""
    if not path_metrics:
        return ""
    
    lines = []
    lines.append("=" * 80)
    lines.append(f"PERFORMANCE RANKING FOR {path} ({START_PHASE_LABEL} -> {PHASE_LABELS.get(path, path)})")
    lines.append("=" * 80)
    
    # 1. Path Quality (DFT Path Avg)
    sorted_path = sorted(path_metrics.items(), key=lambda x: x[1]['avg_path_dft'])
    lines.append("By Path Quality (DFT Avg Energy - Lower is Better):")
    for i, (m, v) in enumerate(sorted_path):
        lines.append(f"  {i+1}. {m:<10} : {v['avg_path_dft']:.5f} eV/at (Barrier: {v['barrier_dft']:.4f})")
    
    # 2. Accuracy (MAE Energy)
    sorted_mae_e = sorted(path_metrics.items(), key=lambda x: x[1]['mae_e'])
    lines.append("\nBy Accuracy (MAE Energy):")
    for i, (m, v) in enumerate(sorted_mae_e):
        lines.append(f"  {i+1}. {m:<10} : {v['mae_e']:.5f} eV/at")
    
    # 3. Stability (Endpoints Force Error)
    sorted_stab = sorted(path_metrics.items(), key=lambda x: x[1]['err_endpts_f'])
    lines.append("\nBy Stability (Endpoint Forces Error):")
    for i, (m, v) in enumerate(sorted_stab):
        lines.append(f"  {i+1}. {m:<10} : {v['err_endpts_f']:.5f} eV/Å")
    
    # 4. Accuracy (MAE Force)
    sorted_mae_f = sorted(path_metrics.items(), key=lambda x: x[1]['mae_f'])
    lines.append("\nBy Accuracy (MAE Force):")
    for i, (m, v) in enumerate(sorted_mae_f):
        lines.append(f"  {i+1}. {m:<10} : {v['mae_f']:.5f} eV/Å")
    
    lines.append("\n" + "=" * 80 + "\n\n")
    return "\n".join(lines)


def apply_phase_ticks(ax, images, end_label):
    """Applies specific Space Group labels to the first and last ticks."""
    tick_labels = [str(x) for x in images]
    tick_labels[0] = START_PHASE_LABEL
    tick_labels[-1] = end_label
    ax.set_xticks(images)
    ax.set_xticklabels(tick_labels, fontsize=FS_TICK)


def plot_comparison(path, mlip, data, limits, output_dir):
    """Create a comparison plot for a single path/MLIP combination."""
    images = data['images']
    ve, vf = data['dft_e'], data['dft_f']
    me, mf = data['mlip_e'], data['mlip_f']
    
    # Identify labels
    target_phase_label = PHASE_LABELS.get(path, path)
    
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
    
    report_text = generate_report_entry(path, mlip, images, me, ve, mf, vf, metrics_dict)
    
    # --- Plotting ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Energy
    ax1.plot(images, ve, marker=MARKERS['DFT'], color=COLORS['DFT'], lw=4, label='DFT', ms=12)
    ax1.plot(images, me, marker=MARKERS.get(mlip, 's'), color=COLORS.get(mlip, 'r'), ls='--', lw=4, label=mlip, ms=12)
    
    # Apply custom ticks
    apply_phase_ticks(ax1, images, target_phase_label)
    
    ax1.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
    # Updated Title: LiF (Start -> End)
    ax1.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label}) - Energy", fontsize=FS_TITLE)
    ax1.tick_params(axis='y', which='major', labelsize=FS_TICK)
    ax1.legend(fontsize=FS_LEGEND, loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
    if path in limits:
        ax1.set_ylim(limits[path]['e_lim'])
    
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
    
    apply_phase_ticks(ax2, images, target_phase_label)
    
    ax2.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
    ax2.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label}) - Forces", fontsize=FS_TITLE)
    ax2.tick_params(axis='y', which='major', labelsize=FS_TICK)
    ax2.legend(fontsize=FS_LEGEND, loc='upper left')
    ax2.grid(True, alpha=0.3)
    if path in limits:
        ax2.set_ylim(limits[path]['f_lim'])
    
    stats_f = (f"Endpts F Err: {err_endpts_f:.4f}\n"
               f"MAE Force:    {mae_f:.4f}")
    
    ax2.text(0.96, 0.96, stats_f, transform=ax2.transAxes, fontsize=FS_TEXT,
             verticalalignment='top', horizontalalignment='right', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.95, edgecolor='#aaaaaa'))
    
    plt.tight_layout()
    
    # Save PNG and SVG
    out_name_base = f"{path}_{mlip}_comparison"
    plt.savefig(os.path.join(output_dir, out_name_base + ".png"), dpi=150)
    plt.savefig(os.path.join(output_dir, out_name_base + ".svg"), format='svg')
    plt.close()
    print(f"Created: {out_name_base}.png & .svg")
    
    return report_text, metrics_dict


def plot_all_mlibs_comparison(path, data_cache, limits, output_dir):
    """Create a combined plot showing all MLIPs vs DFT for a single path."""
    if path not in data_cache:
        return
    
    target_phase_label = PHASE_LABELS.get(path, path)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    
    # Get DFT data (same for all MLIPs)
    first_mlip = list(data_cache[path].keys())[0]
    dft_data = data_cache[path][first_mlip]
    images = dft_data['images']
    ve = dft_data['dft_e']
    vf = dft_data['dft_f']
    
    # Energy plot
    ax1.plot(images, ve, marker=MARKERS['DFT'], color=COLORS['DFT'], lw=4, label='DFT', ms=12)
    for mlip in MLIPS:
        if mlip in data_cache[path]:
            me = data_cache[path][mlip]['mlip_e']
            ax1.plot(images, me, marker=MARKERS.get(mlip, 's'), color=COLORS.get(mlip, 'r'),
                     ls='--', lw=3, label=mlip, ms=10)
    
    apply_phase_ticks(ax1, images, target_phase_label)
    ax1.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
    # New Title
    ax1.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Energy", fontsize=FS_TITLE)
    ax1.tick_params(axis='y', which='major', labelsize=FS_TICK)
    ax1.legend(fontsize=FS_LEGEND - 2, loc='upper left')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
    if path in limits:
        ax1.set_ylim(limits[path]['e_lim'])
    
    # Force plot
    ax2.plot(images, vf, marker=MARKERS['DFT'], color=COLORS['DFT'], lw=4, label='DFT', ms=12)
    for mlip in MLIPS:
        if mlip in data_cache[path]:
            mf = data_cache[path][mlip]['mlip_f']
            ax2.plot(images, mf, marker=MARKERS.get(mlip, 's'), color=COLORS.get(mlip, 'r'),
                     ls='--', lw=3, label=mlip, ms=10)
    
    apply_phase_ticks(ax2, images, target_phase_label)
    ax2.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
    # New Title
    ax2.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Forces", fontsize=FS_TITLE)
    ax2.tick_params(axis='y', which='major', labelsize=FS_TICK)
    ax2.legend(fontsize=FS_LEGEND - 2, loc='upper left')
    ax2.grid(True, alpha=0.3)
    if path in limits:
        ax2.set_ylim(limits[path]['f_lim'])
    
    plt.tight_layout()
    out_name_base = f"{path}_all_MLIPs_comparison"
    plt.savefig(os.path.join(output_dir, out_name_base + ".png"), dpi=150)
    plt.savefig(os.path.join(output_dir, out_name_base + ".svg"), format='svg')
    plt.close()
    print(f"Created: {out_name_base}.png & .svg")


def plot_all_paths_overview(data_cache, output_dir):
    """Create an overview plot showing all paths for each MLIP."""
    for mlip in MLIPS:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
        
        path_colors = {'LiF_Cm': '#1f77b4', 'LiF_P1': '#ff7f0e', 'LiF_Pm': '#2ca02c'}
        
        for path in PATHS:
            if path in data_cache and mlip in data_cache[path]:
                data = data_cache[path][mlip]
                images = data['images']
                ve = data['dft_e']
                vf = data['dft_f']
                me = data['mlip_e']
                mf = data['mlip_f']
                
                # Get clean label
                dest = PHASE_LABELS.get(path, path)
                label_base = f"{START_PHASE_LABEL} $\\rightarrow$ {dest}"

                # Energy - DFT solid, MLIP dashed
                ax1.plot(images, ve, color=path_colors[path], lw=3, label=f"{label_base} DFT")
                ax1.plot(images, me, color=path_colors[path], ls='--', lw=3, label=f"{label_base} {mlip}")
                
                # Force
                ax2.plot(images, vf, color=path_colors[path], lw=3, label=f"{label_base} DFT")
                ax2.plot(images, mf, color=path_colors[path], ls='--', lw=3, label=f"{label_base} {mlip}")
        
        ax1.set_xlabel("Image Index", fontsize=FS_LABEL)
        ax1.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
        ax1.set_title(f"All Paths - Energy ({mlip})", fontsize=FS_TITLE)
        ax1.tick_params(axis='both', which='major', labelsize=FS_TICK)
        ax1.legend(fontsize=FS_LEGEND - 4, loc='upper left', ncol=2)
        ax1.grid(True, alpha=0.3)
        ax1.axhline(0, color='gray', linestyle=':', alpha=0.5)
        
        ax2.set_xlabel("Image Index", fontsize=FS_LABEL)
        ax2.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
        ax2.set_title(f"All Paths - Forces ({mlip})", fontsize=FS_TITLE)
        ax2.tick_params(axis='both', which='major', labelsize=FS_TICK)
        ax2.legend(fontsize=FS_LEGEND - 4, loc='upper left', ncol=2)
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        out_name_base = f"all_paths_{mlip}_overview"
        plt.savefig(os.path.join(output_dir, out_name_base + ".png"), dpi=150)
        plt.savefig(os.path.join(output_dir, out_name_base + ".svg"), format='svg')
        plt.close()
        print(f"Created: {out_name_base}.png & .svg")


def create_summary_csv(data_cache, output_dir):
    """Create a summary CSV with all metrics."""
    csv_path = os.path.join(output_dir, "metrics_summary.csv")
    
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Path', 'MLIP', 'Rxn_E_Err_eV_at', 'Barrier_Err_eV_at',
                         'DFT_Barrier_eV_at', 'MAE_Energy_eV_at', 'MAE_Force_eV_A',
                         'Endpts_F_Err_eV_A', 'Avg_Path_DFT_eV_at'])
        
        for path in PATHS:
            if path in data_cache:
                for mlip in MLIPS:
                    if mlip in data_cache[path]:
                        data = data_cache[path][mlip]
                        ve, vf = data['dft_e'], data['dft_f']
                        me, mf = data['mlip_e'], data['mlip_f']
                        
                        err_rxn = abs(ve[-1] - me[-1])
                        err_barrier = abs(max(ve) - max(me))
                        barrier_dft = max(ve)
                        n = min(len(ve), len(me))
                        mae_e = np.mean(np.abs(ve[1:n] - me[1:n])) if n > 1 else 0.0
                        mae_f = np.mean(np.abs(vf[:n] - mf[:n]))
                        err_endpts_f = abs(vf[0] - mf[0]) + abs(vf[-1] - mf[-1])
                        avg_path_dft = np.mean(ve)
                        
                        writer.writerow([path, mlip, f"{err_rxn:.6f}", f"{err_barrier:.6f}",
                                         f"{barrier_dft:.6f}", f"{mae_e:.6f}", f"{mae_f:.6f}",
                                         f"{err_endpts_f:.6f}", f"{avg_path_dft:.6f}"])
    
    print(f"Created: metrics_summary.csv")


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Loading data...")
    
    # Load DFT data
    dft_file = os.path.join(BASE_DIR, "neb_image_details_DFT.csv")
    dft_data = parse_dft_csv(dft_file)
    
    if not dft_data:
        print("Error: Could not load DFT data")
        return
    
    print(f"Loaded DFT data for paths: {list(dft_data.keys())}")
    
    # Load MLIP data for each path
    data_cache = {}
    
    for path in PATHS:
        data_cache[path] = {}
        path_dir = os.path.join(BASE_DIR, path)
        
        if not os.path.exists(path_dir):
            print(f"Warning: Path directory {path_dir} not found")
            continue
        
        for mlip in MLIPS:
            mlip_file = os.path.join(path_dir, MLIP_FILES[mlip])
            
            # Determine column names based on MLIP
            if mlip == 'MACE-OMAT':
                energy_col = 'mace_omat_energy_ev_per_atom'
                force_col = 'mace_omat_max_force'
            elif mlip == 'MACE-MP':
                energy_col = 'mace_mp_energy_ev_per_atom'
                force_col = 'mace_mp_max_force'
            elif mlip == 'ESEN':
                energy_col = 'esen_energy_ev_per_atom'
                force_col = 'esen_max_force'
            elif mlip == 'UMA':
                energy_col = 'uma_energy_ev_per_atom'
                force_col = 'uma_max_force'
            
            mlip_data = parse_mlip_csv(mlip_file, energy_col, force_col)
            
            if mlip_data and path in dft_data:
                data_cache[path][mlip] = get_aligned_data(dft_data[path], mlip_data)
                print(f"Loaded {mlip} data for {path}")
    
    # Calculate limits
    limits = calculate_limits(data_cache)
    
    # Generate plots and reports
    full_report_content = []
    print("-" * 60)
    
    for path in PATHS:
        path_metrics = {}
        
        # 1. Process Individual MLIPs
        for mlip in MLIPS:
            if mlip in data_cache[path]:
                txt, metrics = plot_comparison(path, mlip, data_cache[path][mlip], limits, OUTPUT_DIR)
                full_report_content.append(txt)
                path_metrics[mlip] = metrics
        
        # 2. Generate combined plot for all MLIPs
        plot_all_mlibs_comparison(path, data_cache, limits, OUTPUT_DIR)
        
        # 3. Generate Ranking for this Path
        ranking_txt = generate_ranking_block(path, path_metrics)
        full_report_content.append(ranking_txt)
    
    # 4. Generate overview plots
    plot_all_paths_overview(data_cache, OUTPUT_DIR)
    
    # 5. Create summary CSV
    create_summary_csv(data_cache, OUTPUT_DIR)
    
    # Save summary report
    summary_path = os.path.join(OUTPUT_DIR, SUMMARY_FILENAME)
    with open(summary_path, 'w') as f:
        f.write("\n".join(full_report_content))
    
    print("-" * 60)
    print(f"Summary saved to: {summary_path}")
    print(f"All plots saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()