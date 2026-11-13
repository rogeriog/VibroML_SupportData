#!/usr/bin/env python3
"""
Merge the corrected 2x2x2 supercell results for Bi2Sn2O7 idx=34 and idx=39
into the five-engine dataset.

The original ip_variants_*.json files used a [1,1,1] supercell fallback for
these structures, producing zero-dispersion (q-independent) phonon spectra.
The corrected runs use a true [2,2,2] supercell.

This script:
1. Loads the original ip_variants_*.json for each engine
2. Loads the corrected ip_variants_*_sc222.json (and idx39-specific files)
3. Replaces the entries for idx=34 and idx=39 with the corrected values
4. Writes merged files: ip_variants_*_merged.json
5. Reports a comparison table

Usage:
  python merge_sc222_results.py
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

ENGINES = ["mace", "macemp", "esen", "uma", "pet"]
CORRECTED_IDXS = [34, 39]


def load_json(path):
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def merge_engine(engine):
    """Merge corrected sc222 results into the original ip_variants file."""
    orig_path = os.path.join(HERE, f"ip_variants_{engine}.json")
    sc222_path = os.path.join(HERE, f"ip_variants_{engine}_sc222.json")
    # eSEN idx=39 may be in a separate file (CPU backup)
    idx39_path = os.path.join(HERE, f"ip_variants_{engine}_sc222_idx39.json")

    orig = load_json(orig_path)
    if orig is None:
        print(f"  SKIP {engine}: original file not found")
        return None

    sc222 = load_json(sc222_path) or []
    idx39_extra = load_json(idx39_path) or []

    # Build a map of corrected results by idx
    corrected = {}
    for r in sc222:
        corrected[r["idx"]] = r
    for r in idx39_extra:
        corrected[r["idx"]] = r

    # Merge: replace entries in orig with corrected ones
    merged = []
    n_replaced = 0
    for r in orig:
        if r["idx"] in corrected:
            merged.append(corrected[r["idx"]])
            n_replaced += 1
        else:
            merged.append(r)

    # Write merged file
    out_path = os.path.join(HERE, f"ip_variants_{engine}_merged.json")
    with open(out_path, "w") as f:
        json.dump(merged, f, indent=2)

    return engine, n_replaced, corrected, orig


def main():
    print("=" * 70)
    print("Merging corrected 2x2x2 supercell results for idx=34 and idx=39")
    print("=" * 70)

    all_results = {}

    for eng in ENGINES:
        result = merge_engine(eng)
        if result is None:
            continue
        engine, n_replaced, corrected, orig = result
        all_results[eng] = (corrected, orig)
        print(f"\n{engine}: replaced {n_replaced}/{len(CORRECTED_IDXS)} entries")

        # Print comparison
        for idx in CORRECTED_IDXS:
            orig_entry = next((r for r in orig if r["idx"] == idx), None)
            new_entry = corrected.get(idx)
            if orig_entry and new_entry:
                old_sf = orig_entry["soft_frac_no_thresh"]
                new_sf = new_entry["soft_frac_no_thresh"]
                old_om = orig_entry["omega_min"]
                new_om = new_entry["omega_min"]
                print(f"  idx={idx}: soft_frac {old_sf:.4f} -> {new_sf:.4f}, "
                      f"omega_min {old_om:.2f} -> {new_om:.2f}")
            elif new_entry:
                print(f"  idx={idx}: NEW (no original), soft_frac={new_entry['soft_frac_no_thresh']:.4f}")
            else:
                print(f"  idx={idx}: NOT CORRECTED (still using original)")

    # Summary table
    print("\n" + "=" * 70)
    print("Summary: corrected soft_frac for idx=34 and idx=39")
    print("=" * 70)
    print(f"{'Engine':<12} {'idx=34 (orig)':>14} {'idx=34 (corr)':>14} "
          f"{'idx=39 (orig)':>14} {'idx=39 (corr)':>14}")
    print("-" * 70)
    for eng in ENGINES:
        if eng not in all_results:
            continue
        corrected, orig = all_results[eng]
        o34 = next((r for r in orig if r["idx"] == 34), None)
        o39 = next((r for r in orig if r["idx"] == 39), None)
        n34 = corrected.get(34)
        n39 = corrected.get(39)
        o34s = f"{o34['soft_frac_no_thresh']:.4f}" if o34 else "—"
        n34s = f"{n34['soft_frac_no_thresh']:.4f}" if n34 else "MISSING"
        o39s = f"{o39['soft_frac_no_thresh']:.4f}" if o39 else "—"
        n39s = f"{n39['soft_frac_no_thresh']:.4f}" if n39 else "MISSING"
        print(f"{eng:<12} {o34s:>14} {n34s:>14} {o39s:>14} {n39s:>14}")

    # Check which engines are missing corrections
    missing = []
    for eng in ENGINES:
        if eng not in all_results:
            continue
        corrected, _ = all_results[eng]
        for idx in CORRECTED_IDXS:
            if idx not in corrected:
                missing.append((eng, idx))

    if missing:
        print(f"\nWARNING: {len(missing)} entries still missing corrections:")
        for eng, idx in missing:
            print(f"  {eng} idx={idx}")

    print("\nMerged files written as ip_variants_*_merged.json")
    print("To use these in the analysis, either:")
    print("  1. Copy merged -> original: cp ip_variants_*_merged.json ip_variants_*.json")
    print("  2. Or point the analysis script at the merged files")


if __name__ == "__main__":
    main()
