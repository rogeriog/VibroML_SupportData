#!/usr/bin/env python3
"""
Generate comparison plots for LiF Phase Transition NEB paths.
FEATURES:
  - Manual Limit Overrides: Hardcoded limits based on user feedback to fix buffers.
  - Broken Axes: Supports split Y-axis for outliers.
  - Scientific Notation Disabled: Forces plain numbers on axes.
  - Limits Output: Prints the exact Y-limits used.
"""

import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

# --- Configuration ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(BASE_DIR, "comparison_plots")

PATHS = ['LiF_Cm', 'LiF_P1', 'LiF_Pm']

PHASE_LABELS = {
    'LiF_Cm': r'$Cm$',
    'LiF_P1': r'$P1$',
    'LiF_Pm': r'$Pm$'
}
START_PHASE_LABEL = r'$P6_3mc$'

MLIPS = ['MACE-OMAT', 'MACE-MP', 'ESEN', 'UMA']

MLIP_FILES = {
    'MACE-OMAT': 'results_mace_omat.csv',
    'MACE-MP': 'results_mace_mp.csv',
    'ESEN': 'results_esen.csv',
    'UMA': 'results_uma.csv'
}

COLORS = {
    'MACE-OMAT': '#d62728', 'MACE-MP': '#ff7f0e',
    'ESEN': '#2ca02c', 'UMA': '#9467bd', 'DFT': '#1f77b4'
}

MARKERS = {
    'MACE-OMAT': 's', 'MACE-MP': 'p',
    'ESEN': '^', 'UMA': 'D', 'DFT': 'o'
}

# --- Font Sizes ---
FS_TITLE = 24
FS_LABEL = 20
FS_TICK = 18
FS_LEGEND = 20
FS_TEXT = 15

# --- MANUAL LIMITS CONFIGURATION ---
# Format: 'Path': { 'metric': { 'type': 'broken'/'standard', 'vals': [ (bot_min, bot_max), (top_min, top_max) ] or (min, max) } }
MANUAL_LIMITS = {
    'LiF_Cm': {
        'energy': {
            'type': 'broken',
            # Bot: widened from -15/12 to -20/15
            # Top: widened from 97.6/97.7 to 96/99 (fixes tight squeeze)
            'vals': [(-20, 15), (96, 99)] 
        },
        'force': {
            'type': 'broken',
            # Bot: widened from -0.5/128 to -10/140
            # Top: widened from 2546.4/2546.5 to 2540/2560 (fixes scientific notation issue)
            'vals': [(-10, 140), (2540, 2560)]
        }
    },
    'LiF_P1': {
        'energy': {
            'type': 'broken',
            # Bot: widened from -4/4 to -8/8
            # Top: widened from 5/36 to 4/45 (more buffer above)
            'vals': [(-8, 8), (8, 45)]
        },
        'force': {
            'type': 'broken',
            # Bot: widened from -0.5/276 to -20/300 (buffer below)
            # Top: widened from 536/2341 to 500/2500 (buffer above/below)
            'vals': [(-20, 300), (400, 2500)]
        }
    },
    'LiF_Pm': {
        'energy': {
            'type': 'broken',
            # Bot: widened from -7.9/5 to -12/6 (buffer below)
            # Top: widened from 8.7/19 to 8/25 (buffer above)
            'vals': [(-12, 6), (8, 25)]
        },
        'force': {
            'type': 'standard',
            # Widened from -0.5/408 to -20/450 (buffer below/above)
            'vals': (-20, 450)
        }
    }
}

# --- Data Parsing Functions ---
def parse_dft_csv(filepath):
    data = {}
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                group = row['Group']
                image_num = int(row['Image'].split('_')[-1])
                if group not in data: data[group] = []
                data[group].append({
                    'image': image_num,
                    'total_energy': float(row['Total_Energy_eV']),
                    'energy_per_atom': float(row['Energy_per_atom_eV']),
                    'max_force': float(row['Max_Force_eV_A'])
                })
        for group in data: data[group].sort(key=lambda x: x['image'])
        return data
    except Exception as e:
        print(f"Error parsing DFT CSV: {e}")
        return None

def parse_mlip_csv(filepath, energy_col, force_col):
    data = []
    if not os.path.exists(filepath): return None
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                image_num = int(row['filename'].split('_')[-1].replace('.cif', ''))
                data.append({
                    'image': image_num,
                    'energy_per_atom': float(row[energy_col]),
                    'max_force': float(row[force_col])
                })
        data.sort(key=lambda x: x['image'])
        return data
    except Exception as e:
        print(f"Error parsing MLIP CSV {filepath}: {e}")
        return None

def get_aligned_data(dft_data, mlip_data):
    dft_e = np.array([d['energy_per_atom'] for d in dft_data])
    dft_f = np.array([d['max_force'] for d in dft_data])
    images = np.array([d['image'] for d in dft_data])
    mlip_e = np.array([d['energy_per_atom'] for d in mlip_data])
    mlip_f = np.array([d['max_force'] for d in mlip_data])
    return {
        'images': images,
        'dft_e': dft_e - dft_e[0],
        'dft_f': dft_f,
        'mlip_e': mlip_e - mlip_e[0],
        'mlip_f': mlip_f
    }

def add_diagonal_break_lines(ax, ax2):
    """Adds the diagonal lines to indicate a broken axis."""
    d = .015 
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((-d, +d), (-d, +d), **kwargs)        
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  

    kwargs.update(transform=ax2.transAxes) 
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  
    ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  

def add_custom_annotations(ax, images, end_label, hide_x_label=False, y_shift=0):
    """Adds Phase labels and Interpolated text."""
    ax.set_xticks(images)
    ax.set_xticklabels(images, fontsize=FS_TICK)
    
    if not hide_x_label:
        trans = ax.get_xaxis_transform()
        y_phase = -0.06 + y_shift
        y_interp = -0.13 + y_shift
        
        ax.text(images[0], y_phase, START_PHASE_LABEL, transform=trans, 
                ha='center', va='top', fontsize=FS_TICK, color='black')
        ax.text(images[-1], y_phase, end_label, transform=trans, 
                ha='center', va='top', fontsize=FS_TICK, color='black')
        
        if len(images) > 2:
            center_x = (images[1] + images[-2]) / 2.0
            ax.text(center_x, y_interp+0.04, "( Interpolated Images )", transform=trans,
                    ha='center', va='top', fontsize=FS_TEXT, color='#444444', style='italic')
        
        ax.set_xlabel("Image Index", fontsize=FS_LABEL, labelpad=35)


# --- Plotting Functions ---

def plot_all_mlibs_comparison(path, data_cache, output_dir):
    """
    Combined plot with MANUAL limits and BROKEN AXIS support.
    """
    if path not in data_cache: return None
    
    target_phase_label = PHASE_LABELS.get(path, path)
    
    first_mlip = list(data_cache[path].keys())[0]
    dft_data = data_cache[path][first_mlip]
    images = dft_data['images']
    
    # Prepare data arrays
    plot_data = [] 
    plot_data.append(('DFT', dft_data['dft_e'], dft_data['dft_f'], COLORS['DFT'], '-'))
    
    for mlip in MLIPS:
        if mlip in data_cache[path]:
            me = data_cache[path][mlip]['mlip_e']
            mf = data_cache[path][mlip]['mlip_f']
            plot_data.append((mlip, me, mf, COLORS.get(mlip, 'r'), '--'))

    # Load Manual Limits
    e_conf = MANUAL_LIMITS.get(path, {}).get('energy', {})
    f_conf = MANUAL_LIMITS.get(path, {}).get('force', {})
    
    is_e_broken = e_conf.get('type') == 'broken'
    is_f_broken = f_conf.get('type') == 'broken'

    fig = plt.figure(figsize=(18, 11)) 
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 3], hspace=0.1, wspace=0.25, top=0.93, bottom=0.18)
    
    # --- Energy Axis Setup ---
    e_lim_info = {}
    if is_e_broken:
        ax_e_top = fig.add_subplot(gs[0, 0])
        ax_e_bot = fig.add_subplot(gs[1, 0], sharex=ax_e_top)
        energy_axes = (ax_e_top, ax_e_bot)
        
        ax_e_top.spines['bottom'].set_visible(False)
        ax_e_bot.spines['top'].set_visible(False)
        ax_e_top.tick_params(labeltop=False, bottom=False)
        ax_e_bot.xaxis.tick_bottom()
        add_diagonal_break_lines(ax_e_top, ax_e_bot)
        
        # Apply Manual Limits
        bot_lims, top_lims = e_conf['vals']
        ax_e_bot.set_ylim(bot_lims)
        ax_e_top.set_ylim(top_lims)
        
        # Disable scientific notation
        ax_e_top.ticklabel_format(useOffset=False, style='plain', axis='y')
        ax_e_bot.ticklabel_format(useOffset=False, style='plain', axis='y')

        e_lim_info = {"type": "broken", "bottom": bot_lims, "top": top_lims}
        
        ax_e_top.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Energy", fontsize=FS_TITLE, pad=45)
        ax_e_bot.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
        
        add_custom_annotations(ax_e_bot, images, target_phase_label)
        plt.setp(ax_e_top.get_xticklabels(), visible=False)
        
    else:
        ax_e_bot = fig.add_subplot(gs[:, 0])
        energy_axes = (None, ax_e_bot)
        ax_e_bot.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Energy", fontsize=FS_TITLE, pad=45)
        ax_e_bot.set_ylabel("Relative Energy (eV/atom)", fontsize=FS_LABEL)
        add_custom_annotations(ax_e_bot, images, target_phase_label)
        
        if 'vals' in e_conf:
            ax_e_bot.set_ylim(e_conf['vals'])
            ax_e_bot.ticklabel_format(useOffset=False, style='plain', axis='y')
            
        e_lim_info = {"type": "standard", "limits": ax_e_bot.get_ylim()}

    # --- Force Axis Setup ---
    f_lim_info = {}
    if is_f_broken:
        ax_f_top = fig.add_subplot(gs[0, 1])
        ax_f_bot = fig.add_subplot(gs[1, 1], sharex=ax_f_top)
        force_axes = (ax_f_top, ax_f_bot)
        
        ax_f_top.spines['bottom'].set_visible(False)
        ax_f_bot.spines['top'].set_visible(False)
        ax_f_top.tick_params(labeltop=False, bottom=False)
        ax_f_bot.xaxis.tick_bottom()
        add_diagonal_break_lines(ax_f_top, ax_f_bot)
        
        # Apply Manual Limits
        bot_lims, top_lims = f_conf['vals']
        ax_f_bot.set_ylim(bot_lims)
        ax_f_top.set_ylim(top_lims)

        # Disable scientific notation explicitly
        ax_f_top.ticklabel_format(useOffset=False, style='plain', axis='y')
        ax_f_bot.ticklabel_format(useOffset=False, style='plain', axis='y')
        
        f_lim_info = {"type": "broken", "bottom": bot_lims, "top": top_lims}
        
        ax_f_top.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Forces", fontsize=FS_TITLE, pad=45)
        ax_f_bot.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
        add_custom_annotations(ax_f_bot, images, target_phase_label)
        plt.setp(ax_f_top.get_xticklabels(), visible=False)
        
    else:
        ax_f_bot = fig.add_subplot(gs[:, 1])
        force_axes = (None, ax_f_bot)
        ax_f_bot.set_title(f"LiF ({START_PHASE_LABEL} $\\rightarrow$ {target_phase_label})\nMLIPs vs DFT Comparison - Forces", fontsize=FS_TITLE, pad=45)
        ax_f_bot.set_ylabel("Max Force (eV/Å)", fontsize=FS_LABEL)
        add_custom_annotations(ax_f_bot, images, target_phase_label, y_shift=0.02)
        
        if 'vals' in f_conf:
            ax_f_bot.set_ylim(f_conf['vals'])
            ax_f_bot.ticklabel_format(useOffset=False, style='plain', axis='y')

        f_lim_info = {"type": "standard", "limits": ax_f_bot.get_ylim()}

    # --- Plotting Loops ---
    # Plot Energies
    for ax in energy_axes:
        if ax is None: continue
        for label, e, f, c, ls in plot_data:
            lw = 4 if label == 'DFT' else 3
            mk = MARKERS.get(label, 'o')
            ax.plot(images, e, marker=mk, color=c, ls=ls, lw=lw, label=label, ms=10)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='gray', linestyle=':', alpha=0.5)
        ax.tick_params(axis='both', labelsize=FS_TICK)

    # Plot Forces
    for ax in force_axes:
        if ax is None: continue
        for label, e, f, c, ls in plot_data:
            lw = 4 if label == 'DFT' else 3
            mk = MARKERS.get(label, 'o')
            ax.plot(images, f, marker=mk, color=c, ls=ls, lw=lw, label=label, ms=10)
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='both', labelsize=FS_TICK)

    # --- LEGEND PLACEMENT (Top Center) ---
    ax_e_legend = energy_axes[0] if energy_axes[0] else energy_axes[1]
    ax_f_legend = force_axes[0] if force_axes[0] else force_axes[1]
    
    ax_e_legend.legend(fontsize=FS_LEGEND, loc='upper left', 
                        bbox_to_anchor=(0.15, 0.90), ncol=1, frameon=False)
    ax_f_legend.legend(fontsize=FS_LEGEND, loc='upper left', 
                        bbox_to_anchor=(0.15, 0.90), ncol=1, frameon=False)

    # Save
    out_name_base = f"{path}_all_MLIPs_comparison"
    plt.savefig(os.path.join(output_dir, out_name_base + ".png"),
            dpi=150, bbox_inches='tight', pad_inches=0.2)

    plt.savefig(os.path.join(output_dir, out_name_base + ".svg"),
            format='svg', bbox_inches='tight', pad_inches=0.2)
    plt.close()
    
    return {"path": path, "energy": e_lim_info, "force": f_lim_info}


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("Loading data...")
    
    dft_file = os.path.join(BASE_DIR, "neb_image_details_DFT.csv")
    dft_data = parse_dft_csv(dft_file)
    if not dft_data: return

    data_cache = {}
    
    for path in PATHS:
        data_cache[path] = {}
        path_dir = os.path.join(BASE_DIR, path)
        if not os.path.exists(path_dir): continue
        
        for mlip in MLIPS:
            mlip_file = os.path.join(path_dir, MLIP_FILES[mlip])
            
            if mlip == 'MACE-OMAT': e_col, f_col = 'mace_omat_energy_ev_per_atom', 'mace_omat_max_force'
            elif mlip == 'MACE-MP': e_col, f_col = 'mace_mp_energy_ev_per_atom', 'mace_mp_max_force'
            elif mlip == 'ESEN': e_col, f_col = 'esen_energy_ev_per_atom', 'esen_max_force'
            elif mlip == 'UMA': e_col, f_col = 'uma_energy_ev_per_atom', 'uma_max_force'
            
            mlip_data = parse_mlip_csv(mlip_file, e_col, f_col)
            if mlip_data and path in dft_data:
                data_cache[path][mlip] = get_aligned_data(dft_data[path], mlip_data)

    print("\n" + "="*80)
    print(f"{'PLOT NAME':<20} | {'METRIC':<10} | {'LIMITS (Y-axis)'}")
    print("-" * 80)

    for path in PATHS:
        result = plot_all_mlibs_comparison(path, data_cache, OUTPUT_DIR)
        if result:
            for metric in ['energy', 'force']:
                info = result[metric]
                label = f"{path} ({metric})"
                if info.get('type') == 'broken':
                    lim_str = f"BOT: {info['bottom']}, TOP: {info['top']} (Broken)"
                else:
                    lim_str = f"{info['limits']}"
                print(f"{label:<33} | {lim_str}")
        
    print("="*80)
    print(f"Plots saved to: {OUTPUT_DIR}\n")

if __name__ == "__main__":
    main()