#!/usr/bin/env python3
"""Methylation tracks of one region from a comet dmr output folder or from per-site BED files.

usage: comet_plot_dmr.py --dmr-dir DMR_OUT {--region chr:start-end | --dmr N [--flank 300]} --out PREFIX
                         [--pool] [--min-depth 1]
       comet_plot_dmr.py --track LABEL=FILE[,FILE...] [--track ...] --region chr:start-end --out PREFIX

--dmr-dir: the --out folder of comet dmr. One track per input file (target files first, then
background files, as listed in 02.run_info/input_files.tsv), labelled with the sample names
without their shared prefix and suffix; the significant DMRs (01.dmr/*.sig.bed) in the region
are shaded. --dmr N plots the N-th significant DMR (row N of *.sig.tsv) with --flank bp on each
side; --pool draws one pooled track per group instead of one per file.
--track: one track per --track (top to bottom), for any per-site BEDs; several comma-separated
FILEs in one track are pooled. Can be added to --dmr-dir tracks.
Each track is the per-CpG methylation profile (n_mod / n_total of CpGs with >= --min-depth reads),
drawn as a line with the area below filled. FILE: bgzipped, tabix-indexed per-site BED (chrom,
start, end, n_mod, n_total), as written by comet assign. Writes PREFIX.pdf and PREFIX.png.
"""
import os
import sys
import csv
import glob
import argparse

import pysam
import matplotlib
matplotlib.use('Agg')
from matplotlib import font_manager
import matplotlib.pyplot as plt

COLOR = '#4C72B0'   # matplotlib / seaborn default blue
DMR_COLOR = '#F2C14E'


def die(msg):
    sys.exit('ERROR: ' + msg)


def use_arial():
    for d in (os.environ.get('COMET_FONT_DIR', ''), os.path.join(sys.prefix, 'fonts'),
              os.path.join(os.environ.get('CONDA_PREFIX', ''), 'fonts')):
        if d and os.path.isdir(d):
            for f in glob.glob(os.path.join(d, '[Aa]rial*.ttf')):
                font_manager.fontManager.addfont(f)
    if 'Arial' in {f.name for f in font_manager.fontManager.ttflist}:
        matplotlib.rcParams['font.family'] = 'Arial'
    else:
        print('WARNING: Arial not found (install mscorefonts or set COMET_FONT_DIR)', file=sys.stderr)
    matplotlib.rcParams.update({'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                                'savefig.dpi': 300, 'savefig.bbox': 'tight'})


def profile(paths, chrom, start, end, min_depth):
    counts = {}
    for path in paths:
        with pysam.TabixFile(path) as tb:
            if chrom not in tb.contigs:
                continue
            for line in tb.fetch(chrom, start, end):
                r = line.split('\t')
                c = counts.setdefault(int(r[1]), [0, 0])
                c[0] += int(r[3])
                c[1] += int(r[4])
    x, y = [], []
    for pos in sorted(counts):
        n_mod, n_tot = counts[pos]
        if n_tot >= min_depth:
            x.append(pos)
            y.append(n_mod / n_tot)
    return x, y


def short_names(names):
    """Sample names without their shared prefix and suffix (brain_deep_AD_1.Microglia.CG -> AD_1)."""
    if len(names) < 2:
        return names
    pre = os.path.commonprefix(names)
    suf = os.path.commonprefix([n[::-1] for n in names])[::-1]
    out = [n[len(pre):len(n) - len(suf)].strip('_.-') for n in names]
    return out if all(out) and len(set(out)) == len(out) else names


def dmr_dir_inputs(d):
    """(sample, group, file) rows of a comet dmr run, target files first."""
    tsv = os.path.join(d, '02.run_info', 'input_files.tsv')
    if not os.path.isfile(tsv):
        die('%s not found (--dmr-dir must be the --out folder of comet dmr)' % tsv)
    with open(tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))
    for r in rows:
        if not os.path.isfile(r['file']):
            die('input BED not found: %s (moved? give the tracks with --track LABEL=FILE)' % r['file'])
    return sorted(rows, key=lambda r: r['group'] != 'target')


def group_names(d):
    """target / background labels given to comet dmr (--target-name, --background-name)."""
    names = {'target': 'target', 'background': 'background'}
    par = os.path.join(d, '02.run_info', 'run_parameters.txt')
    if os.path.isfile(par):
        with open(par) as fh:
            for line in fh:
                f = line.split()
                if len(f) == 2 and f[0] in ('target_name', 'background_name'):
                    names[f[0][:-5]] = f[1]
    return names


def dmr_file(d, ext):
    f = sorted(glob.glob(os.path.join(d, '01.dmr', '*.sig.' + ext)))
    if not f:
        die('no 01.dmr/*.sig.%s in %s' % (ext, d))
    return f[0]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--dmr-dir', help='--out folder of comet dmr')
    p.add_argument('--region', help='chr:start-end (commas allowed)')
    p.add_argument('--dmr', type=int, help='with --dmr-dir: plot the N-th significant DMR (1 = first row of *.sig.tsv)')
    p.add_argument('--flank', type=int, default=300, help='with --dmr: bp added on each side [300]')
    p.add_argument('--pool', action='store_true', help='with --dmr-dir: one pooled track per group')
    p.add_argument('--track', action='append', default=[], help='LABEL=FILE[,FILE...]')
    p.add_argument('--min-depth', type=int, default=1, help='minimum reads per CpG [1]')
    p.add_argument('--out', required=True, help='output prefix (PREFIX.pdf, PREFIX.png)')
    a = p.parse_args()

    if not a.dmr_dir and not a.track:
        die('give --dmr-dir or --track')
    if (a.region is None) == (a.dmr is None):
        die('give either --region or --dmr')
    if a.dmr is not None and not a.dmr_dir:
        die('--dmr needs --dmr-dir')

    if a.region:
        try:
            chrom, rng = a.region.replace(',', '').split(':')
            start, end = (int(v) for v in rng.split('-'))
        except ValueError:
            die('--region must look like chr15:74174700-74175500')
    else:
        with open(dmr_file(a.dmr_dir, 'tsv')) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        if not 1 <= a.dmr <= len(rows):
            die('--dmr must be between 1 and %d (significant DMRs in %s)' % (len(rows), a.dmr_dir))
        r = rows[a.dmr - 1]
        chrom, start, end = r['chrom'], max(0, int(r['start']) - a.flank), int(r['end']) + a.flank
    if end <= start:
        die('region end must be larger than start')

    tracks = []   # (label, [files])
    dmrs = []
    if a.dmr_dir:
        rows = dmr_dir_inputs(a.dmr_dir)
        if a.pool:
            names = group_names(a.dmr_dir)
            for g in dict.fromkeys(r['group'] for r in rows):
                tracks.append((names.get(g, g), [r['file'] for r in rows if r['group'] == g]))
        else:
            for name, r in zip(short_names([r['sample'] for r in rows]), rows):
                tracks.append((name, [r['file']]))
        with open(dmr_file(a.dmr_dir, 'bed')) as fh:
            for line in fh:
                f = line.split('\t')
                if f[0] == chrom and int(f[2]) > start and int(f[1]) < end:
                    dmrs.append((int(f[1]), int(f[2])))
    for t in a.track:
        if '=' not in t:
            die('--track must be LABEL=FILE[,FILE...]: %s' % t)
        label, path = t.split('=', 1)
        tracks.append((label, path.split(',')))

    use_arial()
    n = len(tracks)
    fig, axes = plt.subplots(n, 1, figsize=(3.6, 0.42 * n + 0.4), sharex=True, squeeze=False,
                             gridspec_kw={'hspace': 0.25})
    axes = axes[:, 0]
    for ax, (label, paths) in zip(axes, tracks):
        for s, e in dmrs:
            ax.axvspan(s, e, color=DMR_COLOR, alpha=0.35, lw=0)
        x, y = profile(paths, chrom, start, end, a.min_depth)
        ax.plot(x, y, color=COLOR, lw=0.8)
        ax.fill_between(x, y, 0, color=COLOR, alpha=0.6, lw=0)
        ax.set_ylim(0, 1.2)
        ax.set_yticks([0, 1])
        ax.set_ylabel(label, rotation=0, ha='right', va='center', fontsize=7)
        for side in ('top', 'right'):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis='x', length=0)
        ax.tick_params(axis='y', length=2, labelsize=6)
    last = axes[-1]
    last.set_xlim(start, end)
    ticks = [start, (start + end) // 2, end]
    last.set_xticks(ticks)
    last.set_xticklabels(['{:,}'.format(v) for v in ticks], fontsize=6)
    last.tick_params(axis='x', length=2)
    last.set_xlabel(chrom, fontsize=7)
    d = os.path.dirname(a.out)
    if d:
        os.makedirs(d, exist_ok=True)
    for ext in ('pdf', 'png'):
        fig.savefig('%s.%s' % (a.out, ext))
    print('%s.pdf / .png: %d tracks, %s:%d-%d, %d DMR(s) shaded' % (a.out, n, chrom, start, end, len(dmrs)))


if __name__ == '__main__':
    main()
