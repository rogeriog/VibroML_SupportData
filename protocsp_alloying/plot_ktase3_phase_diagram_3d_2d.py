import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import numpy as np
import warnings

# Suppress runtime warnings for nanmean/nanmax with empty slices
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Read the CSV
df = pd.read_csv('ktase3_phase_diagram_data.csv')

# Get unique polymorphs and separate by dimensionality
polymorphs = df['polymorph'].unique()
top_1d_polys = [p for p in polymorphs if p.startswith('top')]
top_2d_polys = [p for p in polymorphs if p.startswith('2d_top')]

# === Palettes ===
red_orange_shades = [
    '#FF3300', '#FF6347', '#FF8C00', '#FF5733', 
    '#E64A19', '#FF7F50', '#D84315'
]
green_shades = [
    '#008000', '#32CD32', '#00FA9A'
]

poly_colors = {}
idx_1d = 0
idx_2d = 0
for poly in polymorphs:
    if poly == 'base':
        poly_colors[poly] = '#1f77b4'  # Plotly default blue
    elif poly in top_1d_polys:
        poly_colors[poly] = red_orange_shades[idx_1d % len(red_orange_shades)]
        idx_1d += 1
    elif poly in top_2d_polys:
        poly_colors[poly] = green_shades[idx_2d % len(green_shades)]
        idx_2d += 1
    else:
        poly_colors[poly] = '#00CC96'  # Fallback

# === Label Mapping for Legend ===
label_map = {
    'base': '3D perovskite',
    # 1D Chains
    'top1': r'1D-$P2_1/c$ (0.0 meV)',
    'top2': r'1D-$P1$ (+2.8 meV)',
    'top3': r'1D-$P\bar{1}$ (+4.4 meV)',
    'top4': r'1D-$P1$ (+5.5 meV)',
    'top5': r'1D-$P1$ (+6.2 meV)',
    'top6': r'1D-$P2_1/m$ (+7.4 meV)',
    'top7': r'1D-$P1$ (+8.2 meV)',
    # 2D Sheets
    '2d_top1': r'2D-$P1$ (+54.9 meV)',
    '2d_top2': r'2D-$R\bar{3}$ (+62.1 meV)',
    '2d_top3': r'2D-$P1$ (+68.9 meV)',
}

# --- 3D Plotting Function ---
def create_3d_plot(z_column, title, z_title, filename):
    fig = go.Figure()

    # Determine global grid boundaries for unified interpolation
    x_all = df['x'].values
    y_all = df['y'].values
    xi = np.linspace(x_all.min(), x_all.max(), 50)
    yi = np.linspace(y_all.min(), y_all.max(), 50)
    xi, yi = np.meshgrid(xi, yi)

    top_zi_list_1d = []

    # 1. Plot individual surfaces and collect 'top' data
    for polymorph in polymorphs:
        subdf = df[df['polymorph'] == polymorph]
        x = subdf['x'].values
        y = subdf['y'].values
        z = subdf[z_column].values  

        display_name = label_map.get(polymorph, polymorph)

        if len(x) < 4:
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z, 
                mode='markers', 
                name=display_name,
                marker=dict(color=poly_colors[polymorph])
            ))
            continue

        # Interpolate onto the unified grid
        zi = griddata((x, y), z, (xi, yi), method='cubic')
        
        # Store for the volume calculation ONLY if it's a 1D polymorph
        if polymorph in top_1d_polys:
            top_zi_list_1d.append(zi)

        color = poly_colors[polymorph]
        solid_colorscale = [[0, color], [1, color]]

        fig.add_trace(go.Surface(
            x=xi, y=yi, z=zi, 
            name=display_name, 
            colorscale=solid_colorscale, 
            showscale=False,
            showlegend=True 
        ))

    # 2. Calculate and plot the shifted average volume for 1D polymorphs
    if top_zi_list_1d:
        top_zi_array = np.array(top_zi_list_1d)
        z_avg = np.nanmean(top_zi_array, axis=0)
        
        shift_max = np.nanmax(top_zi_array - z_avg)
        shift_min = np.nanmin(top_zi_array - z_avg)
        
        z_upper = z_avg + shift_max
        z_lower = z_avg + shift_min

        light_orange = '#FFD580'
        volume_colorscale = [[0, light_orange], [1, light_orange]]

        fig.add_trace(go.Surface(
            x=xi, y=yi, z=z_upper, 
            name='Envelope Bounds (1D)', 
            colorscale=volume_colorscale, 
            showscale=False,
            opacity=0.3,
            showlegend=True 
        ))

        fig.add_trace(go.Surface(
            x=xi, y=yi, z=z_lower, 
            name='Envelope Bounds (1D) Lower', 
            colorscale=volume_colorscale, 
            showscale=False,
            opacity=0.3,
            showlegend=False 
        ))

    fig.update_layout(
        title=title,
        scene=dict(xaxis_title='x', yaxis_title='y', zaxis_title=z_title),
        width=1000,
        height=800
    )

    fig.write_html(filename)
    print(f"3D plot saved to {filename}")

# --- 2D Heatmap Function ---
def create_2d_heatmap(z_column, title, filename):
    x_all = df['x'].values
    y_all = df['y'].values
    xi = np.linspace(x_all.min(), x_all.max(), 100) # Slightly higher res for 2D
    yi = np.linspace(y_all.min(), y_all.max(), 100)
    xi_mesh, yi_mesh = np.meshgrid(xi, yi)

    # 1. Get the base polymorph grid
    base_df = df[df['polymorph'] == 'base']
    base_zi = griddata((base_df['x'], base_df['y']), base_df[z_column], (xi_mesh, yi_mesh), method='cubic')

    # 2. Get 1D top polymorphs and find the max
    top_zi_list = []
    for poly in top_1d_polys:
        subdf = df[df['polymorph'] == poly]
        if len(subdf['x']) >= 4:
            zi = griddata((subdf['x'], subdf['y']), subdf[z_column], (xi_mesh, yi_mesh), method='cubic')
            top_zi_list.append(zi)
            
    if not top_zi_list:
        print("Not enough 1D 'top' data to generate heatmap.")
        return

    top_zi_array = np.array(top_zi_list)
    top_max_zi = np.nanmax(top_zi_array, axis=0)

    # 3. Calculate difference in meV/atom: (Max(Top) - Base) * 1000
    diff_zi_meV = (top_max_zi - base_zi) * 1000

    # 4. Plot using matplotlib
    plt.figure(figsize=(10, 8))
    
    # Filled contour map
    contour_filled = plt.contourf(xi_mesh, yi_mesh, diff_zi_meV, levels=50, cmap='magma')
    cbar = plt.colorbar(contour_filled)
    cbar.set_label(r'$\Delta E$ (meV/atom)', fontsize=12)
    
    # Add explicit contour lines on top
    contour_lines = plt.contour(xi_mesh, yi_mesh, diff_zi_meV, levels=10, colors='black', alpha=0.4, linewidths=0.8)
    plt.clabel(contour_lines, inline=True, fontsize=9, colors='black', fmt='%1.1f')
    
    # 5. Highlight min and max difference points
    if np.any(~np.isnan(diff_zi_meV)):
        min_idx = np.unravel_index(np.nanargmin(diff_zi_meV), diff_zi_meV.shape)
        max_idx = np.unravel_index(np.nanargmax(diff_zi_meV), diff_zi_meV.shape)
        
        min_x, min_y, min_val = xi_mesh[min_idx], yi_mesh[min_idx], diff_zi_meV[min_idx]
        max_x, max_y, max_val = xi_mesh[max_idx], yi_mesh[max_idx], diff_zi_meV[max_idx]
        
        # Plot markers
        plt.plot(min_x, min_y, marker='o', color='cyan', markersize=10, markeredgecolor='white', label='Max Diff')
        plt.plot(max_x, max_y, marker='^', color='lime', markersize=10, markeredgecolor='white', label='Min Diff')
        
        # Add text annotations (labels swapped as requested)
        plt.text(min_x, min_y - 0.01, f"Max: {min_val:.1f} meV", color='white', 
                 fontsize=11, fontweight='bold', ha='center', va='top', 
                 bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.2'))
        
        plt.text(max_x, max_y + 0.01, f"Min: {max_val:.1f} meV", color='black', 
                 fontsize=11, fontweight='bold', ha='center', va='bottom', 
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'))
                 
        plt.legend(loc='upper right')

    plt.title(title, fontsize=14)
    plt.xlabel('$x$', fontsize=12)
    plt.ylabel('$y$', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"2D Heatmap saved to {filename}")


# --- Execution ---

# 1. 3D Plots
create_3d_plot(
    z_column='formation_energy_per_atom_eV',
    title='KTaSe3 Phase Diagram 3D Surfaces (Formation Energy)',
    z_title='Formation Energy (eV/atom)',
    filename='ktase3_phase_diagram_formation_energy_3d.html'
)

create_3d_plot(
    z_column='total_energy_per_atom_eV',
    title='KTaSe3 Phase Diagram 3D Surfaces (Total Energy)',
    z_title='Total Energy (eV/atom)',
    filename='ktase3_phase_diagram_total_energy_3d.html'
)

# 2. 2D Heatmaps
create_2d_heatmap(
    z_column='formation_energy_per_atom_eV',
    title='Energy Difference: Max(1D) vs 3D perovskite (Formation Energy)',
    filename='ktase3_heatmap_formation_energy.png'
)

create_2d_heatmap(
    z_column='total_energy_per_atom_eV',
    title='Energy Difference: Max(1D) vs 3D perovskite (Total Energy)',
    filename='ktase3_heatmap_total_energy.png'
)