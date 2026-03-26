#!/usr/bin/env python3
"""
Script to compile unique structure counts from all VibroML relaxation runs.

For each run, this script:
1. Finds all summary files across all compound directories
2. Extracts energy per atom values
3. Groups structures by energy (difference < 0.2 meV = 0.0002 eV considered same)
4. Outputs the number of unique structures found by each method

Energy threshold: 0.2 meV = 0.0002 eV
"""

import os
import re
import glob
from pathlib import Path
from collections import defaultdict

# Energy threshold in eV (0.2 meV = 0.0002 eV)
ENERGY_THRESHOLD_EV = 0.0002


def find_all_summary_files(base_dir):
    """
    Find all summary files across all run directories.
    Also identifies runs that need iteration-level processing.
    
    Looks for:
    - overall_relaxation_summary.txt
    - overall_ga_summary.txt  
    - overall_ga_summary_recovered.txt
    - final_phonon_analysis_report.txt
    - iteration-level relaxation_summary.txt (fallback)
    
    Returns: run_data dict, plus a list of runs that need iteration-level processing
    """
    run_data = {}  # run_dir -> {'type': type, 'filepath': filepath}
    runs_needing_iteration = set()  # runs that don't have summary files
    
    # First, find all unique run directories
    all_runs = set()
    for root, dirs, files in os.walk(base_dir):
        parts = [p for p in root.split(os.sep) if p and p != '.']
        for part in parts:
            if 'phonon_output' in part or 'output_' in part:
                # Get the run directory path
                idx = parts.index(part)
                run_dir = os.path.join(*parts[:idx+1])
                all_runs.add(run_dir)
                break
    
    print(f"Found {len(all_runs)} total run directories")
    
    # Now process files to identify which runs have summary files
    for root, dirs, files in os.walk(base_dir):
        for f in files:
            filepath = os.path.join(root, f)
            parts = [p for p in filepath.split(os.sep) if p and p != '.']
            
            # Find run directory
            run_idx = -1
            for i, part in enumerate(parts):
                if 'phonon_output' in part or 'output_' in part:
                    run_idx = i
                    break
            
            if run_idx < 0:
                continue
                
            run_dir = os.path.join(*parts[:run_idx+1])
            
            if f == 'overall_relaxation_summary.txt':
                run_data[run_dir] = {'type': 'overall', 'filepath': filepath}
            elif f == 'overall_ga_summary.txt' or f == 'overall_ga_summary_recovered.txt':
                run_data[run_dir] = {'type': 'ga', 'filepath': filepath}
            elif f == 'final_phonon_analysis_report.txt':
                run_data[run_dir] = {'type': 'opt', 'filepath': filepath}
    
    # Find runs that need iteration-level processing
    for run_dir in all_runs:
        if run_dir not in run_data:
            runs_needing_iteration.add(run_dir)
    
    return run_data, runs_needing_iteration


def parse_ga_recovered_summary(filepath):
    """
    Parse a GA summary file (overall_ga_summary.txt or overall_ga_summary_recovered.txt).
    
    Format:
    Num Atoms    Int. Symbol     Crystal System     Energy (eV/atom)     Iter  Gen   Sample   GA Parameters
    352          Cc              monoclinic         -6.031730               3     3     3        D1:5.30, R21:-0.08, Cell:(0.09,-0.02,-0.08), SC:[2, 2, 2]
    """
    energies = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines
            if not line or line.startswith('---') or line.startswith('Date:') or line.startswith('Total') or line.startswith('All Successful') or line.startswith('Num Atoms'):
                continue
            
            parts = line.split()
            if len(parts) >= 4:
                try:
                    energy = float(parts[3])
                    energies.append(energy)
                except (ValueError, IndexError):
                    continue
    
    return energies


def parse_opt_random_report(filepath):
    """
    Parse a final phonon analysis report file to extract unique structures.
    
    Format:
    Rank  Energy       SpaceGroup   Min Freq     Soft Modes Softest Mode Position
    ----------------------------------------------------------------------------------------------------
    1     -6.031435    P2_12_12_1   -3.0984      2          U @ [0.5, 0.0, 0.5]
    """
    energies = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip header lines and separator lines
            if not line or line.startswith('=') or line.startswith('Rank') or line.startswith('Structure'):
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    # Energy is in the second column
                    energy = float(parts[1])
                    energies.append(energy)
                except (ValueError, IndexError):
                    continue
    
    return energies


def parse_overall_summary(filepath):
    """
    Parse an overall_relaxation_summary.txt file.
    
    Format:
    Num Atoms    Int. Symbol     Crystal System     Energy per Atom (eV/atom) Iter   Sample   Params
    1408         Cc              monoclinic         -6.031678                 0.00 0.00 0.00 0.0 0.0 0.0          3      21       D1:4.000, R21:0.250, Cell:(-0.070, -0.070, -0.070, 0.000, 0.000, 0.000)
    """
    energies = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip header lines and empty lines
            line = line.strip()
            if not line or line.startswith('---') or line.startswith('Analysis') or line.startswith('Num Atoms'):
                continue
            
            # Try to extract energy per atom
            # Format: num_atoms  symbol  crystal_system  energy  ... params
            parts = line.split()
            if len(parts) >= 4:
                try:
                    # Energy is typically the 4th column (index 3)
                    energy = float(parts[3])
                    energies.append(energy)
                except (ValueError, IndexError):
                    continue
    
    return energies


def parse_iteration_summary(filepath):
    """
    Parse an iteration-level relaxation_summary.txt file.
    
    Format varies but typically contains energy information.
    """
    energies = []
    
    try:
        with open(filepath, 'r') as f:
            content = f.read()
            
            # Try various patterns to find energy
            # Pattern 1: "Energy per atom = -X.XXXXXX eV/atom"
            for match in re.finditer(r'Energy per atom\s*=\s*([-+]?\d+\.\d+)', content):
                try:
                    energy = float(match.group(1))
                    if energy < 0:  # Valid energies should be negative
                        energies.append(energy)
                except ValueError:
                    pass
            
            # Pattern 2: Just look for negative float values in reasonable range
            if not energies:
                for match in re.finditer(r'-\d+\.\d+', content):
                    try:
                        energy = float(match.group())
                        if -15 < energy < 0:  # Reasonable energy range for eV/atom
                            energies.append(energy)
                    except ValueError:
                        pass
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
    
    return energies


def group_iteration_files_by_run(iter_summaries):
    """
    Group iteration-level files by their parent run directory.
    Returns dict: run_dir -> list of sample files
    """
    run_groups = defaultdict(list)
    
    for filepath in iter_summaries:
        parts = filepath.split(os.sep)
        
        # Find the run directory (contains 'phonon_output')
        run_idx = -1
        for i, part in enumerate(parts):
            if 'phonon_output' in part or 'output_' in part:
                run_idx = i
                break
        
        if run_idx >= 0:
            run_dir = os.path.join(*parts[:run_idx+1])
        else:
            # Use the first part (compound name) as run_dir
            run_dir = parts[0] if len(parts) > 0 else 'unknown'
        
        run_groups[run_dir].append(filepath)
    
    return run_groups


def count_unique_structures(energies, threshold=ENERGY_THRESHOLD_EV):
    """
    Count unique structures based on energy threshold.
    
    Two structures are considered the same if their energy difference
    is less than the threshold (0.2 meV = 0.0002 eV).
    """
    if not energies:
        return 0, []
    
    # Sort energies
    sorted_energies = sorted(energies)
    
    # Group energies that are within threshold
    unique_groups = []
    current_group = [sorted_energies[0]]
    
    for i in range(1, len(sorted_energies)):
        # Check if this energy is within threshold of the previous one
        if abs(sorted_energies[i] - current_group[-1]) < threshold:
            current_group.append(sorted_energies[i])
        else:
            # Start a new group
            unique_groups.append(current_group)
            current_group = [sorted_energies[i]]
    
    # Don't forget the last group
    unique_groups.append(current_group)
    
    # Each group represents one unique structure
    num_unique = len(unique_groups)
    
    # Get representative energies (e.g., the lowest in each group)
    representative_energies = [min(group) for group in unique_groups]
    
    return num_unique, representative_energies


def get_run_info(filepath):
    """
    Extract compound, method and run directory from the filepath.
    """
    parts = filepath.split(os.sep)
    
    # Filter out empty parts and '.'
    parts = [p for p in parts if p and p != '.']
    
    # Get compound (first directory - should be Bi2Sn2O7, CsPbI3, HfO2, or LiF_simple_cubic)
    compound = parts[0] if parts else 'unknown'
    
    # Find the run directory name
    run_dir = ''
    for part in parts:
        if 'phonon_output' in part or 'output_' in part:
            run_dir = part
            break
    
    if not run_dir:
        run_dir = os.path.basename(os.path.dirname(filepath)) if len(parts) > 1 else 'unknown'
    
    # Determine method from run directory
    method = 'UNKNOWN'
    
    # Check for TRADITIONAL_ALL first (must be before TRADITIONAL check)
    if 'TRADITIONAL_ALL' in run_dir:
        method = 'TRADITIONAL_ALL'
    elif 'phonon_output' in run_dir or 'output_' in run_dir:
        if 'MACE_GA' in run_dir or ('GA' in run_dir and 'MACE' not in run_dir and 'OPT' not in run_dir and 'Traditional' not in run_dir and 'TRADITIONAL' not in run_dir):
            method = 'GA'
        elif 'MACE_TRADITIONAL' in run_dir or 'Traditional' in run_dir:
            method = 'TRADITIONAL'
        elif 'MACE_OPT_RANDOM' in run_dir or 'OPT_RANDOM' in run_dir:
            method = 'OPT_RANDOM'
        elif 'GA5it' in run_dir:
            method = 'GA5it'
    else:
        # Fallback: use full run_dir as method indicator
        if 'GA' in run_dir and 'MACE' not in run_dir:
            method = 'GA'
        elif 'TRADITIONAL' in run_dir:
            method = 'TRADITIONAL'
        elif 'OPT_RANDOM' in run_dir or 'OPT' in run_dir:
            method = 'OPT_RANDOM'
    
    # Add recovered indicator if needed
    if 'recovered' in run_dir.lower():
        method += '_RECOVERED'
    
    return compound, method, run_dir


def count_sample_folders(run_dir):
    """
    Count the total number of sample_* folders in a run directory.
    This represents the total number of structures that were attempted.
    """
    count = 0
    for root, dirs, files in os.walk(run_dir):
        for d in dirs:
            if d.startswith('sample_'):
                count += 1
    return count


def group_iteration_files_by_run(iter_summaries):
    """
    Group iteration-level files by their parent run directory.
    Returns dict: run_dir -> list of sample files
    """
    run_groups = defaultdict(list)
    
    for filepath in iter_summaries:
        parts = [p for p in filepath.split(os.sep) if p and p != '.']
        
        # Find the run directory (contains 'phonon_output')
        run_idx = -1
        for i, part in enumerate(parts):
            if 'phonon_output' in part or 'output_' in part:
                run_idx = i
                break
        
        if run_idx >= 0:
            run_dir = os.path.join(*parts[:run_idx+1])
            run_groups[run_dir].append(filepath)
    
    return run_groups


def process_iteration_run(run_dir, base_dir):
    """
    Process iteration-level relaxation_summary.txt files for a run directory.
    """
    all_energies = []
    
    for root, dirs, files in os.walk(run_dir):
        for f in files:
            if f == 'relaxation_summary.txt':
                filepath = os.path.join(root, f)
                energies = parse_iteration_summary(filepath)
                all_energies.extend(energies)
    
    return all_energies


def main():
    # Search from the current directory (which contains all compound directories)
    base_dir = "."
    output_file = "unique_structures_summary.txt"
    
    # Find all summary files - returns dict keyed by run_dir to avoid duplicates
    run_data, runs_needing_iteration = find_all_summary_files(base_dir)
    
    print(f"Found {len(run_data)} unique run directories with summary files")
    print(f"Found {len(runs_needing_iteration)} run directories needing iteration-level processing\n")
    
    results = []
    
    # Process runs that have summary files
    for run_dir, data in sorted(run_data.items()):
        filepath = data['filepath']
        run_type = data['type']
        
        compound, method, run_name = get_run_info(filepath)
        
        # Parse based on type
        if run_type == 'ga':
            energies = parse_ga_recovered_summary(filepath)
        elif run_type == 'opt':
            energies = parse_opt_random_report(filepath)
        else:
            energies = parse_overall_summary(filepath)
        
        # Count successful relaxations (from summary file)
        num_successful = len(energies)
        
        # Count total attempted (from sample_* folders)
        num_attempted = count_sample_folders(run_dir)
        
        num_unique, representative_energies = count_unique_structures(energies)
        
        results.append({
            'compound': compound,
            'method': method,
            'run_dir': run_name,
            'filepath': filepath,
            'num_successful': num_successful,
            'num_attempted': num_attempted,
            'num_unique': num_unique,
            'energies': representative_energies
        })
        
        print(f"{compound}/{method}: {num_unique} unique (from {num_successful} successful, {num_attempted} attempted)")
        print(f"  Run: {run_name}")
    
    # Process runs that need iteration-level processing
    for run_dir in sorted(runs_needing_iteration):
        compound, method, run_name = get_run_info(run_dir)
        
        energies = process_iteration_run(run_dir, base_dir)
        
        # Count total attempted (from sample_* folders)
        num_attempted = count_sample_folders(run_dir)
        
        if energies:
            num_successful = len(energies)
            num_unique, representative_energies = count_unique_structures(energies)
            
            results.append({
                'compound': compound,
                'method': method,
                'run_dir': run_name,
                'filepath': run_dir,
                'num_successful': num_successful,
                'num_attempted': num_attempted,
                'num_unique': num_unique,
                'energies': representative_energies
            })
            
            print(f"{compound}/{method}: {num_unique} unique (from {num_successful} successful, {num_attempted} attempted)")
            print(f"  Run: {run_name}")
    
    print()
    
    # Write output (simplified - just Unique column)
    with open(output_file, 'w') as f:
        f.write("=" * 90 + "\n")
        f.write("UNIQUE STRUCTURES FOUND BY EACH RUN\n")
        f.write("=" * 90 + "\n\n")
        f.write(f"Energy threshold for uniqueness: {ENERGY_THRESHOLD_EV} eV ({ENERGY_THRESHOLD_EV * 1000:.1f} meV)\n\n")
        
        f.write("-" * 90 + "\n")
        f.write(f"{'Compound':<20} {'Method':<20} {'Run Directory':<45} {'Unique':<8}\n")
        f.write("-" * 90 + "\n")
        
        for r in sorted(results, key=lambda x: (x['compound'], x['method'])):
            f.write(f"{r['compound']:<20} {r['method']:<20} {r['run_dir'][:43]:<45} {r['num_unique']:<8}\n")
        
        f.write("-" * 90 + "\n")
    
    print(f"\nResults written to: {output_file}")


if __name__ == "__main__":
    main()
