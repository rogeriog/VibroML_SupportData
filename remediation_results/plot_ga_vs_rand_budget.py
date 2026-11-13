#!/usr/bin/env python3
"""SI figure: cumulative best energy vs successful relaxations, GA vs RAND.

Data sources, in order of preference per run:
  1. the run's overall summary (``overall_relaxation_summary.txt`` /
     ``overall_ga_summary*.txt`` for GA, ``overall_random_search_summary.txt``
     for RAND) — the complete list of successful relaxations with
     (iter, gen, sample) ordering and per-atom energies;
  2. the per-sample ``relaxation_summary.txt`` files under
     ``main_iter_<I>[_gen_<G>]/sample_<K>/`` when no overall summary exists
     (HfO2 RAND).

Rows whose per-atom energy falls outside the physical window documented
for each compound (steric-clash artifacts, e.g. the 14 LiF RAND rows below
-5.5 eV/atom) are excluded from the curve — the same policy the manuscript
applies to flagged structures elsewhere.

Single-run caveat: each panel shows one GA and one RAND run; outcomes for
the two samplers are mixed and not budget-matched, so the figure does not
establish an advantage of either sampler at fixed compute budget.
"""
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = '/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/02_runs_and_data/VibroML_SupportData/remediation_results'
FIG_DIR = '/gpfs/scratch/acad/htforft/rgouvea/vibroml_runs/03_paper_and_rebuttal/GM_feedback_0911_extracted/VIBROML_PAPER_GMR2/figures'

RUNS = {
    'LiF': {
        'GA': ('LiF_simple_cubic/LiFsimplecubic_GA5it_phonon_output_20250831-152432',
               'overall_relaxation_summary.txt'),
        'RAND': ('LiF_simple_cubic/LiFsimplecubic_MACE_OPT_RANDOM_phonon_output_20251229-142329',
                 'overall_random_search_summary.txt'),
        'e_min': -5.5,            # 14 clash-artifact rows: -6.62 .. -18.57 eV/atom
        'ylim': (-4.85, -4.45),
    },
    'CsPbI$_3$': {
        'GA': ('CsPbI3/CsPbI3_GA5it_phonon_output_20251003-002011',
               'overall_relaxation_summary.txt'),
        'RAND': ('CsPbI3/CsPbI3_MACE_OPT_RANDOM_phonon_output_20251229-001913',
                 'overall_random_search_summary.txt'),
        'e_min': None,
        'ylim': None,
    },
    'HfO$_2$': {
        'GA': ('HfO2/HfO2_Fm3m_GA_phonon_output_20250826-161333',
               'overall_relaxation_summary.txt'),
        'RAND': ('HfO2/HfO2_Fm3m_OPT_RANDOM_phonon_output_20250901-184308',
                 None),                        # no overall summary -> sample dirs
        'e_min': None,
        'ylim': None,
    },
    'Bi$_2$Sn$_2$O$_7$': {
        'GA': ('Bi2Sn2O7/Sn2Bi2O7_MACE_GA_phonon_output_20260107-205131',
               'overall_ga_summary_recovered.txt'),
        'RAND': ('Bi2Sn2O7/Sn2Bi2O7_MACE_OPT_RANDOM_phonon_output_20260123-012721',
                 'overall_random_search_summary.txt'),
        'e_min': None,
        'ylim': None,
    },
}

ROW_RE = re.compile(r'^\s*(\d+)\s+(\S+)\s+(\S+)\s+(-?\d+\.\d+)\s+(.*?)\s*$')
E_RE = re.compile(r'Final Energy per Atom:\s*(-?\d+\.?\d*(?:[Ee][+-]?\d+)?)\s*eV/atom')
N_RE = re.compile(r'Total Number of Atoms:\s*(\d+)')
OK_RE = re.compile(r'Relaxation Status:\s*SUCCESS')
DIR_RE = re.compile(r'main_iter_(\d+)(?:_gen_(\d+))?/sample_(\d+)/relaxation_summary\.txt$')


def parse_overall(path, kind):
    """Parse overall summary rows -> (order_key, energy, natoms).

    GA rows:    natoms sym system energy iter gen sample ...
    RAND rows:  natoms sym system energy [cell params] iter sample ...
    """
    recs = []
    if not os.path.exists(path):
        return recs
    for line in open(path, errors='replace'):
        m = ROW_RE.match(line)
        if not m:
            continue
        natoms = int(m.group(1))
        try:
            energy = float(m.group(4))
        except ValueError:
            continue
        rest = m.group(5).split()
        nums = []
        for tok in rest:
            try:
                nums.append(float(tok))
            except ValueError:
                break
        if kind == 'GA':
            # first three numeric tokens after energy: iter, gen, sample
            if len(nums) < 3:
                continue
            key = (int(nums[0]), int(nums[1]), int(nums[2]))
        else:
            # last two integers before the '{' of random params: iter, sample
            if len(nums) < 2:
                continue
            key = (int(nums[-2]), 0, int(nums[-1]))
        recs.append((key, energy, natoms))
    recs.sort(key=lambda r: r[0])
    return recs


def collect_samples(run_dir):
    """Fallback: per-sample relaxation_summary.txt files."""
    recs = []
    for path in glob.glob(os.path.join(run_dir, 'main_iter_*', '*', 'relaxation_summary.txt')):
        m = DIR_RE.search(path)
        if not m:
            continue
        it, gen, samp = int(m.group(1)), int(m.group(2) or 0), int(m.group(3))
        try:
            txt = open(path, errors='replace').read()
        except OSError:
            continue
        if not OK_RE.search(txt):
            continue
        em, nm = E_RE.search(txt), N_RE.search(txt)
        if not em:
            continue
        recs.append(((it, gen, samp), float(em.group(1)), int(nm.group(1)) if nm else 0))
    recs.sort(key=lambda r: r[0])
    return recs


def get_records(run_dir, summary_file, kind):
    if summary_file:
        recs = parse_overall(os.path.join(run_dir, summary_file), kind)
        if recs:
            return recs
    return collect_samples(run_dir)


def curve(recs, e_min=None):
    """Cumulative-best energy vs successful relaxation index."""
    xs, ys, natoms = [], [], 0
    best = np.inf
    n_plot = 0
    for i, (_, e, n) in enumerate(recs, 1):
        if e_min is not None and e < e_min:
            continue          # steric-clash artifact row
        n_plot += 1
        natoms += n
        if e < best:
            best = e
        xs.append(n_plot)
        ys.append(best)
    return np.array(xs), np.array(ys), len(recs), natoms


fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
axes = axes.ravel()

C_GA, C_RAND = '#0072B2', '#E69F00'   # Wong blue / orange

for ax, (name, spec) in zip(axes, RUNS.items()):
    for method, color in [('GA', C_GA), ('RAND', C_RAND)]:
        rel_dir, summ = spec[method]
        recs = get_records(os.path.join(ROOT, rel_dir), summ, method)
        if not recs:
            print(f'!! no records parsed for {name}/{method}')
            continue
        xs, ys, n_rel, n_at = curve(recs, spec['e_min'])
        n_kept = len(xs)
        note = f' ({n_rel - n_kept} clash-artifact rows excluded)' if n_kept != n_rel else ''
        ax.step(xs, ys, where='post', color=color, linewidth=1.8,
                label=f'{method} ({n_kept} relaxations, {n_at:,} atoms)')
        print(f'{name:14s} {method:4s}: {n_kept:4d}/{n_rel} relaxations, '
              f'{n_at:7,d} atoms, best {ys[-1]:.4f} eV/atom{note}')
    ax.set_title(name, fontsize=14)
    ax.set_xlabel('successful relaxations', fontsize=11)
    ax.set_ylabel('best energy found (eV/atom)', fontsize=11)
    ax.legend(fontsize=9, loc='lower right')
    if spec['ylim']:
        ax.set_ylim(*spec['ylim'])
    ax.tick_params(labelsize=10)

for ext in ('pdf', 'png'):
    out = f'{FIG_DIR}/fig_si_ga_vs_rand_budget.{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print('Saved:', out)
plt.close(fig)
print('Done.')
