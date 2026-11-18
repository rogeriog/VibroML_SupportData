import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import scipy.ndimage as ndi
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.patheffects as path_effects
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings

# === Matplotlib Configuration ===
plt.rcParams['text.usetex'] = False
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['mathtext.fontset'] = 'stix' 

# Suppress runtime warnings for nanmean/nanmax with empty slices
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Read the CSV
df = pd.read_csv('ktase3_phase_diagram_data.csv')

# --- Remove the top7 1D-P1 structure (8.2 meV/at) and 2d_top4 (71.7 meV/at) ---
df = df[~df['polymorph'].isin(['top7', '2d_top4'])]

# Get unique polymorphs
polymorphs = df['polymorph'].unique()
top_1d_polys = [p for p in polymorphs if p.startswith('top') and not p.startswith('2d')]
top_2d_polys = [p for p in polymorphs if p.startswith('2d_top')]
all_polys = ['base'] + top_1d_polys + top_2d_polys

# === Group Definitions for Colors, Legends, and Merging Phase Diagram Regions ===
group_definitions = [
    {'id': 'base', 'polys': ['base'], 'color': '#1f77b4', 'label': '3D perovskite', 'sg': r'3D Perovskite', 'dim': ''},
    {'id': 'top1', 'polys': ['top1'], 'color': '#DC143C', 'label': r'1D-$P2_1/c$ (0.0 meV/at)', 'sg': r'$P2_1/c$', 'dim': '1D'},
    {'id': 'top2_3_4_5', 'polys': ['top2', 'top3', 'top4', 'top5'], 'color': '#FF8C00', 'label': r'1D-$P1$ (+2.8 to +6.2 meV/at)', 'sg': r'$P1$', 'dim': '1D'},
    {'id': 'top6', 'polys': ['top6'], 'color': '#FF69B4', 'label': r'1D-$P2_1/m$ (+7.4 meV/at)', 'sg': r'$P2_1/m$', 'dim': '1D'},
    {'id': '2d_top1', 'polys': ['2d_top1'], 'color': '#32CD32', 'label': r'2D-$P1$ (+54.9 meV/at)', 'sg': r'$P1$', 'dim': '2D'},
    {'id': '2d_top2', 'polys': ['2d_top2'], 'color': '#00FA9A', 'label': r'2D-$R\bar{3}$ (+62.1 meV/at)', 'sg': r'$R\bar{3}$', 'dim': '2D'},
    {'id': '2d_top3', 'polys': ['2d_top3'], 'color': '#2E8B57', 'label': r'2D-$P1$ (+68.9 meV/at)', 'sg': r'$P1$', 'dim': '2D'}
]

# Generate simple lookup dictionary for scatter/surface plotting
poly_colors = {}
for grp in group_definitions:
    for p in grp['polys']:
        poly_colors[p] = grp['color']

def create_5panel_figure(z_column, metric_name, filename_base):
    # Set up the figure for 5 panels in a 3-row layout using GridSpec
    fig = plt.figure(figsize=(24, 27))
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.35, wspace=0.3)
    
    # Stoichiometric title
    title_kt = r'(K$_{1-x}$Ba$_x$)(Ta$_{1-x}$Zr$_x$)(Se$_{1-y}$O$_y$)$_3$'
    fig.suptitle(f'{title_kt}\n{metric_name} Analysis & Phase Diagram', fontsize=28, fontweight='bold', y=0.97)

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
    
    x_center = x_max / 2.0
    y_center = y_max / 2.0

    # High resolution for smooth boundaries in phase diagram
    grid_res = 300
    xi = np.linspace(-x_pad, x_max + x_pad, grid_res)
    yi = np.linspace(-y_pad, y_max + y_pad, grid_res)
    xi_mesh, yi_mesh = np.meshgrid(xi, yi)

    # Pre-compute interpolated grids
    zi_dict = {}
    for polymorph in all_polys:
        subdf = df[df['polymorph'] == polymorph]
        if len(subdf['x']) >= 4:
            zi_dict[polymorph] = griddata((subdf['x'], subdf['y']), subdf[z_column], (xi_mesh, yi_mesh), method='cubic')

    top_zi_array_1d = np.array([zi_dict[p] for p in top_1d_polys if p in zi_dict])
    top_zi_array_2d = np.array([zi_dict[p] for p in top_2d_polys if p in zi_dict])
    all_zi_array = np.array([zi_dict[p] for p in all_polys if p in zi_dict])

    # Keep handles for the global legend
    ordered_handles = []
    for grp in group_definitions:
        patch = mpatches.Patch(color=grp['color'], label=grp['label'])
        ordered_handles.append(patch)

    # ==========================================
    # SUBPLOT 1: 3D Surface Plot (Base vs 1D)
    # ==========================================
    ax1 = fig.add_subplot(gs[0, 0:2], projection='3d')
    for polymorph in ['base'] + top_1d_polys:
        subdf = df[df['polymorph'] == polymorph]
        if len(subdf) == 0: continue
        x, y, z = subdf['x'].values, subdf['y'].values, subdf[z_column].values  
        color = poly_colors[polymorph]
        if polymorph in zi_dict:
            ax1.plot_surface(xi_mesh, yi_mesh, zi_dict[polymorph], color=color, alpha=0.9, rstride=3, cstride=3, antialiased=True)
        else:
            ax1.scatter(x, y, z, color=color, s=70)

    ax1.set_title("Formation Energy\n(3D Perovskite vs 1D Polymorphs)", fontsize=20, pad=15)
    ax1.set_xlabel('x (Ba,Zr fraction)', fontsize=18, labelpad=18)
    ax1.set_ylabel('y (O fraction)', fontsize=18, labelpad=18)
    ax1.set_zlabel('Energy (eV/atom)', fontsize=18, labelpad=18)
    ax1.tick_params(axis='both', which='major', labelsize=14) 
    ax1.view_init(elev=25, azim=-45)

    # ==========================================
    # SUBPLOT 2: 2D Heatmap (Base vs Avg 1D)
    # ==========================================
    ax2 = fig.add_subplot(gs[0, 2:4])
    base_sub = df[df['polymorph'] == 'base'].groupby(['x', 'y'])[z_column].mean()
    top_1d_avg_sub = df[df['polymorph'].isin(top_1d_polys)].groupby(['x', 'y'])[z_column].min()
    
    # Base - Avg(1D) -> Negative means Base is more stable
    exact_diffs_1 = (base_sub - top_1d_avg_sub) * 1000
    diff_zi_1 = (zi_dict['base'] - np.nanmin(top_zi_array_1d, axis=0)) * 1000
    
    vmin_1, vmax_1 = np.nanmin(diff_zi_1), np.nanmax(diff_zi_1)

    contour_1 = ax2.contourf(xi_mesh, yi_mesh, diff_zi_1, levels=50, cmap='magma_r', vmin=vmin_1, vmax=vmax_1)
    ax2.scatter(exact_diffs_1.index.get_level_values('x'), exact_diffs_1.index.get_level_values('y'), 
                c=exact_diffs_1.values, cmap='magma_r', edgecolor='k', s=80, zorder=5, vmin=vmin_1, vmax=vmax_1)

    # Highlight 0.0 line
    ax2.contour(xi_mesh, yi_mesh, diff_zi_1, levels=[0.0], colors='white', linewidths=4.0, linestyles='solid', zorder=8)
    ax2.contour(xi_mesh, yi_mesh, diff_zi_1, levels=[0.0], colors='black', linewidths=1.8, linestyles='dashed', zorder=9)

    min_idx = np.nanargmin(exact_diffs_1.values)
    max_idx = np.nanargmax(exact_diffs_1.values)
    min_x, min_y, min_val = exact_diffs_1.index.get_level_values('x')[min_idx], exact_diffs_1.index.get_level_values('y')[min_idx], exact_diffs_1.values[min_idx]
    max_x, max_y, max_val = exact_diffs_1.index.get_level_values('x')[max_idx], exact_diffs_1.index.get_level_values('y')[max_idx], exact_diffs_1.values[max_idx]

    ax2.plot(min_x, min_y, marker='o', color='cyan', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    ax2.plot(max_x, max_y, marker='^', color='lime', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    
    ha_min = 'left' if min_x <= x_center else 'right'
    va_min = 'bottom' if min_y <= y_center else 'top'
    ha_max = 'left' if max_x <= x_center else 'right'
    va_max = 'bottom' if max_y <= y_center else 'top'
    
    offset = 0.015
    props_min = dict(boxstyle='round,pad=0.3', facecolor='black', alpha=1.0, edgecolor='cyan', linewidth=1.5)
    ax2.text(min_x + (offset if ha_min=='left' else -offset), min_y + (offset if va_min=='bottom' else -offset), 
             f"Min: {min_val:.1f} meV/at", color='white', fontsize=16, fontweight='bold', ha=ha_min, va=va_min, bbox=props_min, zorder=100)
    
    props_max = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=1.0, edgecolor='lime', linewidth=1.5)
    ax2.text(max_x + (offset if ha_max=='left' else -offset), max_y + (offset if va_max=='bottom' else -offset), 
             f"Max: {max_val:.1f} meV/at", color='black', fontsize=16, fontweight='bold', ha=ha_max, va=va_max, bbox=props_max, zorder=100)

    cbar1 = fig.colorbar(contour_1, ax=ax2, fraction=0.046, pad=0.08)
    cbar1.set_label(r'$\Delta E$ (meV/at)', fontsize=18)
    cbar1.ax.tick_params(labelsize=16)
    ax2.clabel(ax2.contour(xi_mesh, yi_mesh, diff_zi_1, levels=10, colors='black', alpha=0.4, linewidths=0.8, zorder=10), inline=True, fontsize=14, colors='black', fmt='%1.1f')

    ax2.set_title("Energy Difference\n3D Perovskite vs lowest 1D polymorph", fontsize=20, pad=15)
    ax2.set_xlabel('x (Ba,Zr fraction)', fontsize=18)
    ax2.set_ylabel('y (O fraction)', fontsize=18)
    ax2.set_xlim(-x_pad, x_max + x_pad); ax2.set_ylim(-y_pad, y_max + y_pad)
    ax2.set_xticks(unique_x); ax2.set_xticklabels(x_labels, fontsize=16)
    ax2.set_yticks(unique_y); ax2.set_yticklabels(y_labels, fontsize=16)

    # ==========================================
    # SUBPLOT 3: 3D Surface Plot (1D vs 2D)
    # ==========================================
    ax3 = fig.add_subplot(gs[1, 0:2], projection='3d')
    for polymorph in top_1d_polys + top_2d_polys:
        subdf = df[df['polymorph'] == polymorph]
        if len(subdf) == 0: continue
        x, y, z = subdf['x'].values, subdf['y'].values, subdf[z_column].values  
        color = poly_colors[polymorph]
        if polymorph in zi_dict:
            ax3.plot_surface(xi_mesh, yi_mesh, zi_dict[polymorph], color=color, alpha=0.9, rstride=3, cstride=3, antialiased=True)
        else:
            ax3.scatter(x, y, z, color=color, s=70)

    ax3.set_title("Formation Energy\n(1D vs 2D Polymorphs)", fontsize=20, pad=15)
    ax3.set_xlabel('x (Ba,Zr fraction)', fontsize=18, labelpad=18)
    ax3.set_ylabel('y (O fraction)', fontsize=18, labelpad=18)
    ax3.set_zlabel('Energy (eV/atom)', fontsize=18, labelpad=18)
    ax3.tick_params(axis='both', which='major', labelsize=14) 
    ax3.view_init(elev=25, azim=-45)

    # ==========================================
    # SUBPLOT 4: 2D Heatmap (Avg 2D vs Avg 1D) -> REVERSED!
    # ==========================================
    ax4 = fig.add_subplot(gs[1, 2:4])
    top_2d_avg_sub = df[df['polymorph'].isin(top_2d_polys)].groupby(['x', 'y'])[z_column].min()
    
    # REVERSED: Avg(2D) - Avg(1D) -> Negative means 2D is more stable
    # Magma_r makes negative values light (yellow/white), indicating 2D stability
    exact_diffs_2 = (top_2d_avg_sub - top_1d_avg_sub) * 1000
    diff_zi_2 = (np.nanmin(top_zi_array_2d, axis=0) - np.nanmin(top_zi_array_1d, axis=0)) * 1000
    
    vmin_2, vmax_2 = np.nanmin(diff_zi_2), np.nanmax(diff_zi_2)

    contour_2 = ax4.contourf(xi_mesh, yi_mesh, diff_zi_2, levels=50, cmap='magma_r', vmin=vmin_2, vmax=vmax_2)
    ax4.scatter(exact_diffs_2.index.get_level_values('x'), exact_diffs_2.index.get_level_values('y'), 
                c=exact_diffs_2.values, cmap='magma_r', edgecolor='k', s=80, zorder=5, vmin=vmin_2, vmax=vmax_2)

    # Highlight 0.0 line
    ax4.contour(xi_mesh, yi_mesh, diff_zi_2, levels=[0.0], colors='white', linewidths=4.0, linestyles='solid', zorder=8)
    ax4.contour(xi_mesh, yi_mesh, diff_zi_2, levels=[0.0], colors='black', linewidths=1.8, linestyles='dashed', zorder=9)

    min_idx2 = np.nanargmin(exact_diffs_2.values)
    max_idx2 = np.nanargmax(exact_diffs_2.values)
    min_x2, min_y2, min_val2 = exact_diffs_2.index.get_level_values('x')[min_idx2], exact_diffs_2.index.get_level_values('y')[min_idx2], exact_diffs_2.values[min_idx2]
    max_x2, max_y2, max_val2 = exact_diffs_2.index.get_level_values('x')[max_idx2], exact_diffs_2.index.get_level_values('y')[max_idx2], exact_diffs_2.values[max_idx2]

    ax4.plot(min_x2, min_y2, marker='o', color='cyan', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    ax4.plot(max_x2, max_y2, marker='^', color='lime', markersize=14, markeredgecolor='white', markeredgewidth=2, zorder=6)
    
    ha_min2 = 'left' if min_x2 <= x_center else 'right'
    va_min2 = 'bottom' if min_y2 <= y_center else 'top'
    ha_max2 = 'left' if max_x2 <= x_center else 'right'
    va_max2 = 'bottom' if max_y2 <= y_center else 'top'

    ax4.text(min_x2 + (offset if ha_min2=='left' else -offset), min_y2 + (offset if va_min2=='bottom' else -offset), 
             f"Min: {min_val2:.1f} meV/at", color='white', fontsize=16, fontweight='bold', ha=ha_min2, va=va_min2, bbox=props_min, zorder=100)
    
    ax4.text(max_x2 + (offset if ha_max2=='left' else -offset), max_y2 + (offset if va_max2=='bottom' else -offset), 
             f"Max: {max_val2:.1f} meV/at", color='black', fontsize=16, fontweight='bold', ha=ha_max2, va=va_max2, bbox=props_max, zorder=100)

    cbar2 = fig.colorbar(contour_2, ax=ax4, fraction=0.046, pad=0.08)
    cbar2.set_label(r'$\Delta E$ (meV/at)', fontsize=18)
    cbar2.ax.tick_params(labelsize=16)
    ax4.clabel(ax4.contour(xi_mesh, yi_mesh, diff_zi_2, levels=10, colors='black', alpha=0.4, linewidths=0.8, zorder=10), inline=True, fontsize=14, colors='black', fmt='%1.1f')

    ax4.set_title("Energy Difference\nlowest 2D polymorph vs lowest 1D polymorph", fontsize=20, pad=15)
    ax4.set_xlabel('x (Ba,Zr fraction)', fontsize=18)
    ax4.set_ylabel('y (O fraction)', fontsize=18)
    ax4.set_xlim(-x_pad, x_max + x_pad); ax4.set_ylim(-y_pad, y_max + y_pad)
    ax4.set_xticks(unique_x); ax4.set_xticklabels(x_labels, fontsize=16)
    ax4.set_yticks(unique_y); ax4.set_yticklabels(y_labels, fontsize=16)

    # ==========================================
    # SUBPLOT 5: 2D Phase Diagram (Bottom Center)
    # ==========================================
    ax5 = fig.add_subplot(gs[2, 1:3])
    
    # 1. Compute Raw Min Phase Index
    min_phase_idx_raw = np.full(xi_mesh.shape, np.nan)
    valid_mask = ~np.isnan(all_zi_array).all(axis=0)
    min_phase_idx_raw[valid_mask] = np.nanargmin(all_zi_array[:, valid_mask], axis=0)

    # 2. Map raw index to our Grouped Index to naturally merge adjacent areas with same label/color
    min_phase_idx_grouped = np.full(xi_mesh.shape, np.nan)
    for group_idx, grp in enumerate(group_definitions):
        for poly in grp['polys']:
            if poly in all_polys:
                raw_idx = all_polys.index(poly)
                min_phase_idx_grouped[min_phase_idx_raw == raw_idx] = group_idx

    cmap_colors = [grp['color'] for grp in group_definitions]
    levels = np.arange(len(group_definitions) + 1) - 0.5
    
    ax5.contourf(xi_mesh, yi_mesh, min_phase_idx_grouped, levels=levels, colors=cmap_colors, alpha=0.8)
    ax5.contour(xi_mesh, yi_mesh, min_phase_idx_grouped, levels=levels, colors='k', linewidths=1.5, alpha=0.7)
    
    ax5.scatter(x_all, y_all, c='black', s=15, alpha=0.5, zorder=5)

    # Render Text dynamically for each Group
    for group_idx, grp in enumerate(group_definitions):
        mask = (min_phase_idx_grouped == group_idx)
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
                # This calculates the distance from every pixel in the blob to the nearest edge
                dist = ndi.distance_transform_edt(largest_mask)
                
                # Get the row/col index of the maximum distance (the deepest point inside the blob)
                max_dist_idx = np.unravel_index(np.argmax(dist), dist.shape)
                
                # Extract the exact physical coordinates for that point
                cx = xi_mesh[max_dist_idx]
                cy = yi_mesh[max_dist_idx]
                
                # Print Main SG Name stacking upwards
                txt_sg = ax5.text(cx, cy, grp['sg'], fontsize=24, fontweight='bold', 
                               ha='center', va='bottom', color='black', zorder=10)
                txt_sg.set_path_effects([path_effects.withStroke(linewidth=4, foreground='white')])
                
                # Print Dimensionality stacking downwards
                if grp['dim']:
                    txt_dim = ax5.text(cx, cy, f"({grp['dim']})", fontsize=16, fontweight='bold', 
                                   ha='center', va='top', color='black', zorder=10)
                    txt_dim.set_path_effects([path_effects.withStroke(linewidth=4, foreground='white')])

    ax5.set_title("Phase Diagram\n(Most Stable Polymorph)", fontsize=20, pad=15)
    ax5.set_xlabel('x (Ba,Zr fraction)', fontsize=18)
    ax5.set_ylabel('y (O fraction)', fontsize=18)
    ax5.set_xlim(-x_pad, x_max + x_pad); ax5.set_ylim(-y_pad, y_max + y_pad)
    ax5.set_xticks(unique_x); ax5.set_xticklabels(x_labels, fontsize=16)
    ax5.set_yticks(unique_y); ax5.set_yticklabels(y_labels, fontsize=16)
    
    # Global legend attached to the left of the Phase Diagram
    ax5.legend(handles=ordered_handles, loc='center right', bbox_to_anchor=(-0.15, 0.5), fontsize=15, ncol=1, title="Phases", title_fontsize=17)

    # Final Layout and Save
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    
    jpg_filename = f'{filename_base}.jpg'
    plt.savefig(jpg_filename, dpi=300, bbox_inches='tight')
    print(f"Saved: {jpg_filename}")
    
    svg_filename = f'{filename_base}.svg'
    plt.savefig(svg_filename, format='svg', bbox_inches='tight')
    print(f"Saved: {svg_filename}")
    
    plt.close()

# --- Execution ---
create_5panel_figure(
    z_column='formation_energy_per_atom_eV',
    metric_name='Formation Energy',
    filename_base='ktase3_combined_formation_energy_5panel'
)