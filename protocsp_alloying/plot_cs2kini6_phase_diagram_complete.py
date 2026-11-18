import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import scipy.ndimage as ndi
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
import warnings
import math

# === Matplotlib Configuration ===
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['mathtext.fontset'] = 'stix' 

# Suppress runtime warnings for nanmean/nanmax with empty slices
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Read the CSV
df = pd.read_csv('cs2kini6_phase_diagram_data.csv')

# Get unique polymorphs and sort them properly (base, then top1 -> top7)
polymorphs = df['polymorph'].unique()
# Sort the 'top' polymorphs explicitly by their numeric suffix
top_polys = sorted([p for p in polymorphs if p.startswith('top')], key=lambda x: int(x.replace('top', '')))
all_polys = ['base'] + top_polys

# === Palettes ===
colors = [
    '#DC143C', # Crimson Red
    '#FFD700', # Yellow/Gold
    '#FF69B4', # Hot Pink
    '#FF4500', # Orange Red
    '#FF7F50', # Coral
    '#C71585', # Medium Violet Red
    '#FF6347'  # Tomato
]

poly_colors = {'base': '#1f77b4'} # Keep 3D perovskite blue
for i, poly in enumerate(top_polys):
    poly_colors[poly] = colors[i % len(colors)]

# === Label Mapping for Legend & Text ===
label_map = {
    'base': r'3D perovskite ($Fm\bar{3}m$)',
    'top1': r'$I\bar{4}2m$ (0.00 meV/at)',
    'top2': r'$P\bar{3}1m$ (+2.27 meV/at)',
    'top3': r'$C2/c$ (+3.25 meV/at)',
    'top4': r'$P\bar{3}$ (+4.07 meV/at)',
    'top5': r'$P\bar{1}$ (+4.58 meV/at)',
    'top6': r'$Fmm2$ (+5.72 meV/at)',
    'top7': r'$P2_1/c$ (+8.30 meV/at)'
}

# Short space group names for the phase diagram regions
sg_map = {
    'base': r'$Fm\bar{3}m$ (3D)',
    'top1': r'$I\bar{4}2m$',
    'top2': r'$P\bar{3}1m$',
    'top3': r'$C2/c$',
    'top4': r'$P\bar{3}$',
    'top5': r'$P\bar{1}$',
    'top6': r'$Fmm2$',
    'top7': r'$P2_1/c$'
}

def create_3panel_figure(z_column, metric_name, filename_base):
    # Set up the figure for 3 panels in a 2-row layout using GridSpec
    fig = plt.figure(figsize=(24, 18))
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.4, wspace=0.3)
    
    # Stoichiometric title
    title_kt = r'Cs$_2$(K$_{1-x}$Na$_x$)In(I$_{1-y}$Br$_y$)$_6$'
    fig.suptitle(f'{title_kt}\n{metric_name} Analysis & Phase Diagram', fontsize=26, fontweight='bold', y=0.98)

    # Grid preparation with margins
    x_all = df['x'].values
    y_all = df['y'].values
    
    unique_x = np.sort(np.unique(np.round(x_all, 2)))
    unique_y = np.sort(np.unique(np.round(y_all, 2)))
    x_labels = [f"{val:.2f}" for val in unique_x]
    y_labels = [f"{val:.2f}" for val in unique_y]

    x_max, y_max = x_all.max(), y_all.max()
    x_pad = x_max * 0.03
    y_pad = y_max * 0.03
    
    # Center points for dynamic label alignment
    x_center = x_max / 2.0
    y_center = y_max / 2.0

    # Use a HIGH resolution for smooth boundaries in the phase diagram
    grid_res = 300
    xi = np.linspace(-x_pad, x_max + x_pad, grid_res)
    yi = np.linspace(-y_pad, y_max + y_pad, grid_res)
    xi_mesh, yi_mesh = np.meshgrid(xi, yi)

    # Pre-compute interpolated grids for all polymorphs
    zi_dict = {}
    for polymorph in polymorphs:
        subdf = df[df['polymorph'] == polymorph]
        if len(subdf['x']) >= 4:
            zi_dict[polymorph] = griddata((subdf['x'], subdf['y']), subdf[z_column], (xi_mesh, yi_mesh), method='cubic')

    # Compute array for calculations
    top_zi_array = np.array([zi_dict[p] for p in top_polys if p in zi_dict])
    all_zi_array = np.array([zi_dict[p] for p in all_polys if p in zi_dict])

    # ==========================================
    # SUBPLOT 1: 3D Surface Plot (Top Left)
    # ==========================================
    ax1 = fig.add_subplot(gs[0, 0:2], projection='3d')
    
    # Keep the raw ordered list for the legend in panel 3
    ordered_handles = []

    for polymorph in all_polys:
        subdf = df[df['polymorph'] == polymorph]
        x, y, z = subdf['x'].values, subdf['y'].values, subdf[z_column].values  
        color = poly_colors[polymorph]
        display_label = label_map.get(polymorph, polymorph)
        patch = mpatches.Patch(color=color, label=display_label)
        
        ordered_handles.append(patch)

        if polymorph in zi_dict:
            ax1.plot_surface(xi_mesh, yi_mesh, zi_dict[polymorph], color=color, alpha=0.9, rstride=3, cstride=3, antialiased=True)
        else:
            ax1.scatter(x, y, z, color=color, s=70, label=display_label)

    ax1.set_title("Formation Energy\n(3D Perovskite vs 7 Polymorphs)", fontsize=20, pad=15)
    ax1.set_xlabel('x (Na fraction)', fontsize=18, labelpad=18)
    ax1.set_ylabel('y (Br fraction)', fontsize=18, labelpad=18)
    ax1.set_zlabel('Energy (eV/atom)', fontsize=18, labelpad=18)
    ax1.tick_params(axis='both', which='major', labelsize=14) 
    ax1.view_init(elev=25, azim=-45)

    # ==========================================
    # SUBPLOT 2: 2D Heatmap (Top Right)
    # ==========================================
    ax2 = fig.add_subplot(gs[0, 2:4])
    base_sub = df[df['polymorph'] == 'base'].groupby(['x', 'y'])[z_column].mean()
    top_avg_sub = df[df['polymorph'].isin(top_polys)].groupby(['x', 'y'])[z_column].min()
    
    # Inverted: 3D Perovskite - Avg(7 Polymorphs)
    exact_diffs = (base_sub - top_avg_sub) * 1000
    diff_zi = (zi_dict['base'] - np.nanmin(top_zi_array, axis=0)) * 1000
    
    vmin, vmax = np.nanmin(diff_zi), np.nanmax(diff_zi)

    contour = ax2.contourf(xi_mesh, yi_mesh, diff_zi, levels=50, cmap='magma_r', vmin=vmin, vmax=vmax)
    ax2.scatter(exact_diffs.index.get_level_values('x'), exact_diffs.index.get_level_values('y'), 
                c=exact_diffs.values, cmap='magma_r', edgecolor='k', s=80, zorder=5, vmin=vmin, vmax=vmax)

    # --- NEW: Enhanced 0.0 Contour Line Highlight ---
    # First, draw a thick solid white line as a base background
    ax2.contour(xi_mesh, yi_mesh, diff_zi, levels=[0.0], colors='white', linewidths=4.0, linestyles='solid', zorder=8)
    # Then, draw a thinner dashed black line directly over it for high contrast
    ax2.contour(xi_mesh, yi_mesh, diff_zi, levels=[0.0], colors='black', linewidths=1.8, linestyles='dashed', zorder=9)
    # -----------------------------------------------

    min_idx = np.nanargmin(exact_diffs.values)
    max_idx = np.nanargmax(exact_diffs.values)
    
    min_x, min_y, min_val = exact_diffs.index.get_level_values('x')[min_idx], exact_diffs.index.get_level_values('y')[min_idx], exact_diffs.values[min_idx]
    max_x, max_y, max_val = exact_diffs.index.get_level_values('x')[max_idx], exact_diffs.index.get_level_values('y')[max_idx], exact_diffs.values[max_idx]

    ax2.plot(min_x, min_y, marker='o', color='cyan', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    ax2.plot(max_x, max_y, marker='^', color='lime', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    
    ha_min = 'left' if min_x <= x_center else 'right'
    va_min = 'bottom' if min_y <= y_center else 'top'
    ha_max = 'left' if max_x <= x_center else 'right'
    va_max = 'bottom' if max_y <= y_center else 'top'
    
    offset = 0.015
    off_x_min = offset if ha_min == 'left' else -offset
    off_y_min = offset if va_min == 'bottom' else -offset
    off_x_max = offset if ha_max == 'left' else -offset
    off_y_max = offset if va_max == 'bottom' else -offset

    props_min = dict(boxstyle='round,pad=0.3', facecolor='black', alpha=1.0, edgecolor='cyan', linewidth=1.5)
    ax2.text(min_x + off_x_min, min_y + off_y_min, f"Min: {min_val:.1f} meV/at", color='white', fontsize=16, fontweight='bold', ha=ha_min, va=va_min, bbox=props_min, zorder=12)
    
    props_max = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=1.0, edgecolor='lime', linewidth=1.5)
    ax2.text(max_x + off_x_max, max_y + off_y_max, f"Max: {max_val:.1f} meV/at", color='black', fontsize=16, fontweight='bold', ha=ha_max, va=va_max, bbox=props_max, zorder=12)

    cbar = fig.colorbar(contour, ax=ax2, fraction=0.046, pad=0.08)
    cbar.set_label(r'$\Delta E$ (meV/at)', fontsize=18)
    cbar.ax.tick_params(labelsize=16)
    #clabel contours left on for generic levels, labels added on top of highlight
    ax2.clabel(ax2.contour(xi_mesh, yi_mesh, diff_zi, levels=10, colors='black', alpha=0.4, linewidths=0.8, zorder=10), inline=True, fontsize=14, colors='black', fmt='%1.1f')

    ax2.set_title("Energy Difference\n3D perovskite vs lowest competing polymorph", fontsize=20, pad=15)
    ax2.set_xlabel('x (Na fraction)', fontsize=18)
    ax2.set_ylabel('y (Br fraction)', fontsize=18)
    ax2.set_xlim(-x_pad, x_max + x_pad); ax2.set_ylim(-y_pad, y_max + y_pad)
    ax2.set_xticks(unique_x); ax2.set_xticklabels(x_labels, fontsize=16)
    ax2.set_yticks(unique_y); ax2.set_yticklabels(y_labels, fontsize=16)

    # ==========================================
    # SUBPLOT 3: 2D Phase Diagram (Bottom Center)
    # ==========================================
    ax3 = fig.add_subplot(gs[1, 1:3])
    
    min_phase_idx = np.full(xi_mesh.shape, np.nan)
    valid_mask = ~np.isnan(all_zi_array).all(axis=0)
    min_phase_idx[valid_mask] = np.nanargmin(all_zi_array[:, valid_mask], axis=0)

    cmap_colors = [poly_colors[p] for p in all_polys]
    levels = np.arange(len(all_polys) + 1) - 0.5
    
    ax3.contourf(xi_mesh, yi_mesh, min_phase_idx, levels=levels, colors=cmap_colors, alpha=0.8)
    ax3.contour(xi_mesh, yi_mesh, min_phase_idx, levels=levels, colors='k', linewidths=1.5, alpha=0.7)
    
    ax3.scatter(x_all, y_all, c='black', s=15, alpha=0.5, zorder=5)

    for i, polymorph in enumerate(all_polys):
        mask = (min_phase_idx == i)
        if np.any(mask):
            # 1. Find connected components (blobs) for this phase
            labeled_array, num_features = ndi.label(mask)
            
            if num_features > 0:
                # Count the size of each blob (ignore index 0, which is the background)
                sizes = np.bincount(labeled_array.ravel())[1:]
                
                # Identify the label of the largest blob
                largest_blob_label = sizes.argmax() + 1
                
                # Create a new mask containing ONLY the largest blob
                largest_mask = (labeled_array == largest_blob_label)
                
                # 2. Find the "thickest" part of the blob to guarantee the label fits inside
                dist = ndi.distance_transform_edt(largest_mask)
                
                # Get the row/col index of the maximum distance (deepest point)
                max_dist_idx = np.unravel_index(np.argmax(dist), dist.shape)
                
                # Extract the exact physical coordinates for that point
                cx = xi_mesh[max_dist_idx]
                cy = yi_mesh[max_dist_idx]
                
                txt = ax3.text(cx, cy, sg_map[polymorph], fontsize=28, fontweight='bold', 
                               ha='center', va='center', color='black', zorder=10)
                txt.set_path_effects([path_effects.withStroke(linewidth=5, foreground='white')])
                
    ax3.set_title("Phase Diagram\n(Most Stable Polymorph)", fontsize=20, pad=15)
    ax3.set_xlabel('x (Na fraction)', fontsize=18)
    ax3.set_ylabel('y (Br fraction)', fontsize=18)
    ax3.set_xlim(-x_pad, x_max + x_pad); ax3.set_ylim(-y_pad, y_max + y_pad)
    ax3.set_xticks(unique_x); ax3.set_xticklabels(x_labels, fontsize=16)
    ax3.set_yticks(unique_y); ax3.set_yticklabels(y_labels, fontsize=16)
    
    # Legend correctly uses the raw, energy-sorted handles array
    ax3.legend(handles=ordered_handles, loc='center right', bbox_to_anchor=(-0.15, 0.5), fontsize=15, ncol=1, title="Phases", title_fontsize=17)

    # Final Layout and Save
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    
    jpg_filename = f'{filename_base}.jpg'
    plt.savefig(jpg_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {jpg_filename}")
    
    svg_filename = f'{filename_base}.svg'
    plt.savefig(svg_filename, format='svg', bbox_inches='tight')
    print(f"Saved: {svg_filename}")
    
    plt.close()

# --- Execution ---
create_3panel_figure(
    z_column='formation_energy_per_atom_eV',
    metric_name='Formation Energy',
    filename_base='cs2kini6_combined_formation_energy_3panel'
)