#!/usr/bin/env python3
"""
Script to compile timing information from VibroML calculations across different compounds and modes.
"""

import os
import re
import json
import gzip
from datetime import datetime
from pathlib import Path
import csv

# Base directory
BASE_DIR = Path("/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/examples_VibroML/VibroML_SupportData")

# Define compounds and their modes
COMPOUNDS = {
    "Bi2Sn2O7": {
        "GA": ["Sn2Bi2O7_MACE_GA_phonon_output_20260107-205131"],
        "OPT_RANDOM": ["Sn2Bi2O7_MACE_OPT_RANDOM_phonon_output_20260123-012721"],
        "TRADITIONAL": ["Sn2Bi2O7_MACE_TRADITIONAL_phonon_output_20251229-145533"],
        "TRADITIONAL_ALL": ["Sn2Bi2O7_MACE_TRADITIONAL_ALL_phonon_output_20260123-012716"]
    },
    "CsPbI3": {
        "GA": ["CsPbI3_MACE_GA_phonon_output_20251229-135624"],
        "OPT_RANDOM": ["CsPbI3_MACE_OPT_RANDOM_phonon_output_20251229-001913"],
        "TRADITIONAL": ["CsPbI3_TRADITIONAL_phonon_output_20251002-224829"],
        "TRADITIONAL_ALL": ["CsPbI3_TRADITIONAL_ALL_phonon_output_20251003-170449"]
    },
    "HfO2": {
        "GA": ["HfO2_Fm3m_GA_phonon_output_20250826-161333"],
        "OPT_RANDOM": ["HfO2_Fm3m_OPT_RANDOM_phonon_output_20250901-184308"],
        "TRADITIONAL": ["HfO2_Fm3m_MACE_TRADITIONAL_phonon_output_20260108-162322"],
        "TRADITIONAL_ALL": ["HfO2_Fm3m_TRADITIONAL_ALL_phonon_output_20250903-162016"]
    },
    "LiF_simple_cubic": {
        "GA": ["LiFsimplecubic_GA5it_phonon_output_20250831-152432"],
        "OPT_RANDOM": ["LiFsimplecubic_MACE_OPT_RANDOM_phonon_output_20251229-142329"],
        "TRADITIONAL": ["LiFsimplecubic_Traditional_run_phonon_output_20250706-173614"],
        "TRADITIONAL_ALL": ["LiFsimplecubic_TRADITIONAL_ALL_phonon_output_20250904-003516"]
    }
}

def parse_analysis_date(content):
    """Extract Analysis Date from file content."""
    match = re.search(r"Analysis Date:\s+(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})", content)
    if match:
        return datetime.strptime(match.group(1), "%Y-%m-%d %H:%M:%S")
    return None

def find_first_timing_ga(run_dir):
    """Find timing for first relaxation in GA mode."""
    paths_to_try = [
        run_dir / "main_iter_1_gen_1" / "sample_1" / "relaxation_summary.txt",
        run_dir / "main_iter_1" / "sample_1" / "relaxation_summary.txt",
    ]
    for path in paths_to_try:
        if path.exists():
            content = path.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(path)
    return None, None

def find_last_timing_ga(run_dir):
    """Find timing for final generation in GA mode."""
    main_iter_dirs = []
    for d in run_dir.iterdir():
        if d.is_dir() and d.name.startswith("main_iter_") and "_gen_" in d.name:
            if (d / "relaxation_summary_generation.txt").exists():
                main_iter_dirs.append(d)
    
    if not main_iter_dirs:
        return None, None
    
    def sort_key(x):
        parts = x.name.replace("main_iter_", "").split("_gen_")
        if len(parts) < 2:
            return (-1, -1)
        try:
            iter_num = int(parts[0])
            gen_num = int(parts[1])
            return (iter_num, gen_num)
        except ValueError:
            return (-1, -1)
    
    main_iter_dirs.sort(key=sort_key)
    last_iter_dir = main_iter_dirs[-1]
    
    gen_file = last_iter_dir / "relaxation_summary_generation.txt"
    if gen_file.exists():
        content = gen_file.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(gen_file)
    
    return None, None

def find_first_timing_opt_random(run_dir):
    """Find timing for first relaxation in OPT_RANDOM mode."""
    first_iter = run_dir / "main_iter_1"
    if first_iter.exists():
        iter_summary = first_iter / "iteration_1_summary.txt"
        if iter_summary.exists():
            content = iter_summary.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(iter_summary)
        
        sample_file = first_iter / "sample_1" / "relaxation_summary.txt"
        if sample_file.exists():
            content = sample_file.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(sample_file)
    
    return None, None

def find_last_timing_opt_random(run_dir):
    """Find timing for final iteration in OPT_RANDOM mode."""
    iter_dirs = []
    for d in run_dir.iterdir():
        if d.is_dir() and d.name.startswith("main_iter_"):
            try:
                iter_num = int(d.name.replace("main_iter_", ""))
                iter_dirs.append((iter_num, d))
            except:
                pass
    
    if not iter_dirs:
        return None, None
    
    iter_dirs.sort(key=lambda x: x[0])
    last_iter_dir = iter_dirs[-1][1]
    iter_num = iter_dirs[-1][0]
    
    iter_summary = last_iter_dir / f"iteration_{iter_num}_summary.txt"
    if iter_summary.exists():
        content = iter_summary.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(iter_summary)
    
    sample_file = last_iter_dir / "sample_1" / "relaxation_summary.txt"
    if sample_file.exists():
        content = sample_file.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(sample_file)
    
    return None, None

def find_first_timing_traditional(run_dir):
    """Find timing for first relaxation in TRADITIONAL mode."""
    initial_dir = run_dir / "initial_relaxation_for_single_run"
    if initial_dir.exists():
        relax_files = list(initial_dir.glob("relaxation_summary*.txt"))
        if relax_files:
            content = relax_files[0].read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(relax_files[0])
    
    first_iter = run_dir / "iter_1"
    if first_iter.exists():
        relax_iter = first_iter / "relaxation_summary_iter.txt"
        if relax_iter.exists():
            content = relax_iter.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(relax_iter)
        
        sample_file = first_iter / "sample_1" / "relaxation_summary.txt"
        if sample_file.exists():
            content = sample_file.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(sample_file)
    
    return None, None

def find_last_timing_traditional(run_dir):
    """Find timing for final iteration in TRADITIONAL mode."""
    iter_dirs = []
    for d in run_dir.iterdir():
        if d.is_dir() and d.name.startswith("iter_"):
            rest = d.name[5:]
            if rest.isdigit():
                iter_dirs.append(d)
    
    if not iter_dirs:
        overall_file = run_dir / "overall_relaxation_summary.txt"
        if overall_file.exists():
            content = overall_file.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(overall_file)
        return None, None
    
    iter_dirs.sort(key=lambda x: int(x.name[5:]))
    last_iter_dir = iter_dirs[-1]
    
    relax_file = last_iter_dir / "relaxation_summary_iter.txt"
    if relax_file.exists():
        content = relax_file.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(relax_file)
    
    relax_files = list(last_iter_dir.glob("*/relaxation_summary.txt"))
    if relax_files:
        content = relax_files[0].read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(relax_files[0])
    
    return None, None

def find_first_timing_traditional_all(run_dir):
    """Find timing for first iteration in TRADITIONAL_ALL mode."""
    paths_to_try = [
        run_dir / "iter_1" / "iteration_1_summary.txt",
    ]
    for path in paths_to_try:
        if path.exists():
            content = path.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(path)
    return None, None

def find_last_timing_traditional_all(run_dir):
    """Find timing for final iteration in TRADITIONAL_ALL mode."""
    iter_dirs = []
    for d in run_dir.iterdir():
        if d.is_dir() and d.name.startswith("iter_"):
            rest = d.name[5:]
            if rest.isdigit():
                iter_dirs.append(d)
    
    if not iter_dirs:
        summary_file = run_dir / "overall_traditional_all_summary.txt"
        if summary_file.exists():
            content = summary_file.read_text()
            date = parse_analysis_date(content)
            if date:
                return date, str(summary_file)
        return None, None
    
    iter_dirs.sort(key=lambda x: int(x.name[5:]))
    last_iter_dir = iter_dirs[-1]
    iter_num = int(last_iter_dir.name[5:])
    
    summary_file = last_iter_dir / f"iteration_{iter_num}_summary.txt"
    if summary_file.exists():
        content = summary_file.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(summary_file)
    
    overall_file = run_dir / "overall_traditional_all_summary.txt"
    if overall_file.exists():
        content = overall_file.read_text()
        date = parse_analysis_date(content)
        if date:
            return date, str(overall_file)
    
    return None, None

def get_checkpoint_timing(run_dir, duration_hours):
    """Calculate actual running time by subtracting only gaps >= 6 hours."""
    checkpoint_dir = run_dir / "checkpoints"
    if not checkpoint_dir.exists():
        return None, []
    
    timestamps = []
    for cf in checkpoint_dir.iterdir():
        if cf.name.startswith("checkpoint_") and cf.name.endswith(".json.gz"):
            if cf.name in ["checkpoint_metadata.json.gz", "checkpoint_latest.json.gz"]:
                continue
            parts = cf.name.replace("checkpoint_", "").replace(".json.gz", "").split("_")
            if len(parts) >= 2:
                try:
                    dt = datetime.strptime(parts[0] + parts[1], "%Y%m%d%H%M%S")
                    timestamps.append(dt.timestamp())
                except:
                    pass
    
    if len(timestamps) < 2:
        return None, []
    
    timestamps.sort()
    gap_threshold = 21600  # 6 hours
    gaps = []
    total_gaps_seconds = 0
    
    for i in range(1, len(timestamps)):
        gap = timestamps[i] - timestamps[i-1]
        if gap >= gap_threshold:
            gaps.append({
                "start": datetime.fromtimestamp(timestamps[i-1]).strftime("%Y-%m-%d %H:%M"), 
                "end": datetime.fromtimestamp(timestamps[i]).strftime("%Y-%m-%d %H:%M"), 
                "gap_hours": round(gap/3600, 2)
            })
            total_gaps_seconds += gap
    
    # Calculate actual running time: Wall Duration - Sum of detected gaps
    actual_running_time_hours = round(duration_hours - (total_gaps_seconds / 3600), 2)
    
    return {
        "actual_running_time_hours": actual_running_time_hours,
        "gaps": gaps
    }, timestamps

def get_timing_functions(mode):
    """Return the appropriate timing extraction functions for a mode."""
    if mode == "GA":
        return find_first_timing_ga, find_last_timing_ga
    elif mode == "OPT_RANDOM":
        return find_first_timing_opt_random, find_last_timing_opt_random
    elif mode == "TRADITIONAL":
        return find_first_timing_traditional, find_last_timing_traditional
    elif mode == "TRADITIONAL_ALL":
        return find_first_timing_traditional_all, find_last_timing_traditional_all
    else:
        return None, None

def main():
    results = []
    
    for compound, modes in COMPOUNDS.items():
        for mode, runs in modes.items():
            for run_name in runs:
                run_dir = BASE_DIR / compound / run_name
                
                if not run_dir.exists():
                    print(f"Warning: Run directory not found: {run_dir}")
                    continue
                
                first_fn, last_fn = get_timing_functions(mode)
                if first_fn is None:
                    print(f"Unknown mode: {mode}")
                    continue
                
                first_date, first_path = first_fn(run_dir)
                last_date, last_path = last_fn(run_dir)
                
                if first_date and last_date:
                    duration = (last_date - first_date).total_seconds() / 3600  # Wall hours
                    
                    # Get checkpoint timing
                    checkpoint_info, _ = get_checkpoint_timing(run_dir, duration)
                    
                    actual_running_hours = duration
                    gap_info = ""
                    
                    if checkpoint_info:
                        # Only use the adjusted actual time if gaps were actually found
                        if checkpoint_info["gaps"]:
                            actual_running_hours = checkpoint_info["actual_running_time_hours"]
                            gap_details = [f"{g['gap_hours']}h from {g['start']} to {g['end']}" for g in checkpoint_info["gaps"]]
                            gap_info = f" (gaps: {', '.join(gap_details)})"
                    
                    results.append({
                        "compound": compound,
                        "mode": mode,
                        "run_dir": f"./{compound}/{run_name}",
                        "start_date": first_date.strftime("%Y-%m-%d %H:%M:%S"),
                        "end_date": last_date.strftime("%Y-%m-%d %H:%M:%S"),
                        "duration_hours": round(duration, 2),
                        "actual_running_hours": round(actual_running_hours, 2),
                        "start_path": first_path,
                        "end_path": last_path
                    })
                    
                    print(f"{compound} {mode}: {first_date.strftime('%Y-%m-%d %H:%M')} to {last_date.strftime('%Y-%m-%d %H:%M')} = {duration:.2f}h (wall) / {actual_running_hours:.2f}h (actual){gap_info}")
                else:
                    print(f"Warning: Could not find timing for {compound} {mode} ({run_name})")
    
    # Write to CSV
    output_file = BASE_DIR / "vibroml_timing_table.csv"
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["compound", "mode", "run_dir", "start_date", "end_date", "duration_hours", "actual_running_hours", "start_path", "end_path"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\nResults written to {output_file}")

if __name__ == "__main__":
    main()