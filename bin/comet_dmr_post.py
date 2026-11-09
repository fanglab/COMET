#!/usr/bin/env python3
"""Standardise DMR-tool output and apply the comet dmr filters.

standardize TOOL RAW OUT_BED [--alpha A]
    raw metilene output -> sorted BED-like table:
    chrom, start, end, tool_n_cpg, mean_target, mean_background, diff, p_value, q_value,
    statistic, tool_significant   (diff = target - background, methylation as fraction)

finalize STD_BED MAPPED SAMPLE_SHEET OUT_PREFIX --min-cpg N --min-depth D --min-diff X --names T,B
         [--dmr-level group|individual]
    group:      |mean_target - mean_background| >= X
    individual: additionally every target sample differs from every background sample by >= X
                in the direction of the group difference (sample level = summed n_mod / n_total)
    MAPPED: STD_BED rows followed by the number of CpGs covered in any sample and, per sample in
    sample-sheet order, the summed n_mod and n_total over the region; each sample's methylation
    level over the region (summed n_mod / n_total) is written as mean_<sample>.
    Writes OUT_PREFIX.all.tsv, OUT_PREFIX.sig.tsv and OUT_PREFIX.sig.bed.
"""
import sys
import argparse

STD_COLS = ['chrom', 'start', 'end', 'tool_n_cpg', 'mean_target', 'mean_background', 'diff',
            'p_value', 'q_value', 'statistic', 'tool_significant']


def standardize(a):
    rows = []
    with open(a.raw) as f:
        lines = [l.rstrip('\n').split('\t') for l in f if l.strip() and not l.startswith('##')]
    if a.tool == 'metilene':
        # chr start end q mean_diff nCpG MWU_p KS_p mean_a mean_b  (a = target, b = background)
        for r in lines:
            if r[0] in ('chr', 'chrom'):
                continue
            q = float(r[3])
            rows.append([r[0], r[1], r[2], r[5], r[8], r[9], r[4], r[6], r[3], r[7],
                         int(q < a.alpha)])
    rows.sort(key=lambda x: (x[0], int(x[1]), int(x[2])))
    with open(a.out_bed, 'w') as out:
        for r in rows:
            out.write('\t'.join(str(v) for v in r) + '\n')


def finalize(a):
    samples = [l.rstrip('\n').split('\t') for l in open(a.sample_sheet) if l.strip()]
    is_target = [s[1] == 'target' for s in samples]
    tname, bname = a.names.split(',')
    hdr = ['chrom', 'start', 'end', 'n_cpg', 'n_cpg_any_sample', 'mean_' + tname, 'mean_' + bname,
           'diff_%s_minus_%s' % (tname, bname), 'p_value', 'q_value', 'statistic',
           'tool_significant', 'min_mean_depth', 'min_individual_diff', 'pass_dmr_diff',
           'pass_individual_diff', 'pass_min_cpg', 'pass_min_depth', 'significant']
    hdr += ['mean_' + s[0] for s in samples]
    hdr += ['meth_' + s[0] for s in samples] + ['total_' + s[0] for s in samples]
    n_all = n_sig = 0
    with open(a.mapped) as f, open(a.out_prefix + '.all.tsv', 'w') as fa, \
            open(a.out_prefix + '.sig.tsv', 'w') as fs, open(a.out_prefix + '.sig.bed', 'w') as fb:
        fa.write('\t'.join(hdr) + '\n')
        fs.write('\t'.join(hdr) + '\n')
        k = len(STD_COLS)
        for line in f:
            r = line.rstrip('\n').split('\t')
            std, rest = r[:k], r[k:]
            # n_cpg: the CpG count reported by the DMR tool (as in the manuscript's depth filter);
            # n_cpg_any_sample: CpGs in the region covered in at least one sample
            n_any = int(rest[0]) if rest[0] not in ('.', '') else 0
            n_cpg = int(std[3]) if std[3] not in ('NA', '') else n_any
            vals = [0 if v in ('.', '') else int(float(v)) for v in rest[1:]]
            meth, total = vals[0::2], vals[1::2]
            depth = [t / n_cpg if n_cpg else 0.0 for t in total]
            diff = float(std[6])
            pass_diff = int(abs(diff) >= a.min_diff - 1e-12)
            # smallest target-vs-background gap over all sample pairs, signed by the group direction
            lv = [m / t if t else None for m, t in zip(meth, total)]
            lt = [v for v, it in zip(lv, is_target) if it]
            lb = [v for v, it in zip(lv, is_target) if not it]
            if None in lt or None in lb or not lt or not lb:
                gap = None
            else:
                gap = (min(lt) - max(lb)) if diff >= 0 else (min(lb) - max(lt))
            pass_ind = int(gap is not None and gap >= a.min_diff - 1e-12)
            pass_cpg = int(n_cpg >= a.min_cpg)
            pass_depth = int(bool(depth) and min(depth) >= a.min_depth)
            pass_level = pass_diff and (pass_ind or a.dmr_level == 'group')
            significant = int(std[10] == '1' and pass_level and pass_cpg and pass_depth)
            out = std[:3] + [str(n_cpg), str(n_any)] + std[4:] + ['%.2f' % (min(depth) if depth else 0),
                                                       'NA' if gap is None else '%.6f' % gap,
                                                       str(pass_diff), str(pass_ind), str(pass_cpg),
                                                       str(pass_depth), str(significant)]
            out += ['NA' if v is None else '%.6f' % v for v in lv]
            out += [str(v) for v in meth] + [str(v) for v in total]
            fa.write('\t'.join(out) + '\n')
            n_all += 1
            if significant:
                fs.write('\t'.join(out) + '\n')
                fb.write('\t'.join(std[:3] + ['%s:%s' % (tname, std[6])]) + '\n')
                n_sig += 1
    print('DMRs reported by tool: %d | significant after dmr-diff / min-cpg / min-depth filters: %d'
          % (n_all, n_sig))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('standardize')
    s.add_argument('tool', choices=['metilene'])
    s.add_argument('raw')
    s.add_argument('out_bed')
    s.add_argument('--alpha', type=float, default=0.05)
    f = sub.add_parser('finalize')
    f.add_argument('std_bed')
    f.add_argument('mapped')
    f.add_argument('sample_sheet')
    f.add_argument('out_prefix')
    f.add_argument('--min-cpg', type=int, required=True)
    f.add_argument('--min-depth', type=float, required=True)
    f.add_argument('--min-diff', type=float, required=True)
    f.add_argument('--dmr-level', choices=['group', 'individual'], default='group')
    f.add_argument('--names', required=True)
    a = p.parse_args()
    standardize(a) if a.cmd == 'standardize' else finalize(a)


if __name__ == '__main__':
    main()
