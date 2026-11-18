import os
import json
import csv
import re

def extract_data():
    base_path = "lemat_unique_dataset/ProtoCSP"
    data = []

    # Folders to process
    folders = ["cs2kini6_phase_diagram"] + [f"cs2kini6_phase_diagram_top{i}" for i in range(1, 8)]

    for folder in folders:
        folder_path = os.path.join(base_path, folder)
        if not os.path.exists(folder_path):
            continue

        # Find all *_mlip_energies.json files recursively
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                if file.endswith("_mlip_energies.json"):
                    file_path = os.path.join(root, file)
                    with open(file_path, 'r') as f:
                        json_data = json.load(f)

                    # Find rank 1
                    rank_1 = next((item for item in json_data if item.get("rank") == 1), None)
                    if rank_1:
                        # Extract x and y from path
                        path_parts = root.split(os.sep)
                        x_match = re.search(r'x_(\d+p\d+)', root)
                        y_match = re.search(r'y_(\d+p\d+)', root)
                        if x_match and y_match:
                            x_str = x_match.group(1).replace('p', '.')
                            y_str = y_match.group(1).replace('p', '.')
                            x = float(x_str)
                            y = float(y_str)
                        else:
                            continue  # Skip if no x/y

                        # Determine polymorph
                        if folder == "cs2kini6_phase_diagram":
                            polymorph = "base"
                        else:
                            polymorph = folder.split("_")[-1]  # e.g., top1

                        row = {
                            "base_compound": "Cs2KInI6",
                            "formula": rank_1["formula"],
                            "x": x,
                            "y": y,
                            "polymorph": polymorph,
                            "formation_energy_per_atom_eV": rank_1["formation_energy_per_atom_eV"],
                            "total_energy_per_atom_eV": rank_1["total_energy_per_atom_eV"]
                        }
                        data.append(row)

    # Sort data by polymorph, then x, then y
    data.sort(key=lambda row: (row['polymorph'], row['x'], row['y']))

    # Write to CSV
    with open("cs2kini6_phase_diagram_data.csv", 'w', newline='') as csvfile:
        fieldnames = ["base_compound", "formula", "x", "y", "polymorph", "formation_energy_per_atom_eV", "total_energy_per_atom_eV"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

if __name__ == "__main__":
    extract_data()