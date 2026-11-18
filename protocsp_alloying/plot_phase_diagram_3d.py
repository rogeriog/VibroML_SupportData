import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pc
from scipy.interpolate import griddata
import numpy as np

# Read the CSV
df = pd.read_csv('phase_diagram_data.csv')

# Get unique polymorphs and assign a color to each
polymorphs = df['polymorph'].unique()
colors = pc.qualitative.Plotly  # Default Plotly color palette
poly_colors = {poly: colors[i % len(colors)] for i, poly in enumerate(polymorphs)}

def create_3d_plot(z_column, title, z_title, filename):
    # Create figure
    fig = go.Figure()

    # For each polymorph, create a surface
    for polymorph in polymorphs:
        subdf = df[df['polymorph'] == polymorph]
        x = subdf['x'].values
        y = subdf['y'].values
        z = subdf[z_column].values  

        if len(x) < 4:  # Not enough points for surface, skip or plot as scatter
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z, 
                mode='markers', 
                name=polymorph,
                marker=dict(color=poly_colors[polymorph])
            ))
            continue

        # Create grid for interpolation
        xi = np.linspace(x.min(), x.max(), 50)
        yi = np.linspace(y.min(), y.max(), 50)
        xi, yi = np.meshgrid(xi, yi)
        zi = griddata((x, y), z, (xi, yi), method='cubic')

        # Define a solid colorscale for the current polymorph
        color = poly_colors[polymorph]
        solid_colorscale = [[0, color], [1, color]]

        # Add surface
        fig.add_trace(go.Surface(
            x=xi, y=yi, z=zi, 
            name=polymorph, 
            colorscale=solid_colorscale, 
            showscale=False,
            showlegend=True  # Enables the legend for the surface trace
        ))

    # Update layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='x',
            yaxis_title='y',
            zaxis_title=z_title
        ),
        width=1000,
        height=800
    )

    # Save to HTML
    fig.write_html(filename)
    print(f"3D plot saved to {filename}")

# 1. Plot for Formation Energy
create_3d_plot(
    z_column='formation_energy_per_atom_eV',
    title='Phase Diagram 3D Surfaces (Formation Energy)',
    z_title='Formation Energy (eV/atom)',
    filename='phase_diagram_formation_energy_3d.html'
)

# 2. Plot for Total Energy
create_3d_plot(
    z_column='total_energy_per_atom_eV',
    title='Phase Diagram 3D Surfaces (Total Energy)',
    z_title='Total Energy (eV/atom)',
    filename='phase_diagram_total_energy_3d.html'
)