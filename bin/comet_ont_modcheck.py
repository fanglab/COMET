#!/usr/bin/env python3
"""Check which cytosine modifications an ONT modBAM carries before running COMET.

Samples primary autosomal reads with MM/ML tags and reports:
  5mC present, 5hmC present, and whether calls are CpG-only (e.g. dorado 5mCG / 5mCG_5hmCG
  models) or all-context (5mC / 5mC_5hmC models, calls also at CH cytosines).
Exits non-zero when the requested outputs (--mods) are not supported by the data.
"""
import re
import sys
import argparse

import pysam

AUTOSOME = re.compile(r'^chr([1-9]|1[0-9]|2[0-2])$')
CH_FRACTION_ALL_CONTEXT = 0.05   # CpG-only models give ~0 calls at CH cytosines


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bam')
    p.add_argument('--mods', default='5mCG', help='comma list of 5mCG,5hmCG,5mCH,5hmCH')
    p.add_argument('--n-reads', type=int, default=2000)
    a = p.parse_args()

    pysam.set_verbosity(0)   # silence htslib warnings about individual malformed MM tags
    n_reads = 0
    calls = {'m': [0, 0], 'h': [0, 0]}      # code -> [CG calls, CH calls]
    with pysam.AlignmentFile(a.bam, 'rb') as bam:
        for chrom in [c for c in bam.references if AUTOSOME.match(c)]:
            for read in bam.fetch(chrom):
                if read.flag & 0x904 or not read.has_tag('MM'):
                    continue
                seq = read.get_forward_sequence()
                mods = read.modified_bases_forward or {}
                for (base, _strand, code), pairs in mods.items():
                    if base != 'C' or code not in calls:
                        continue
                    for qpos, _ml in pairs:
                        cg = qpos + 1 < len(seq) and seq[qpos + 1] == 'G'
                        calls[code][0 if cg else 1] += 1
                n_reads += 1
                if n_reads >= a.n_reads:
                    break
            if n_reads >= a.n_reads:
                break

    has_m = sum(calls['m']) > 0
    has_h = sum(calls['h']) > 0
    tot_c = sum(calls['m']) + sum(calls['h'])
    ch_frac = (calls['m'][1] + calls['h'][1]) / tot_c if tot_c else 0.0
    all_context = ch_frac > CH_FRACTION_ALL_CONTEXT
    kind = ('5mC' if all_context else '5mCG') + (('_5hmC' if all_context else '_5hmCG') if has_h else '')
    sys.stderr.write('modBAM check (%d reads): 5mC=%s 5hmC=%s context=%s (CH fraction %.3f) -> %s\n' % (
        n_reads, 'yes' if has_m else 'no', 'yes' if has_h else 'no',
        'all' if all_context else 'CpG-only', ch_frac, kind if has_m else 'unsupported'))

    errors = []
    if n_reads == 0:
        errors.append('no primary autosomal reads with MM/ML tags were found')
    elif not has_m:
        errors.append('no 5mC (C+m) calls found; COMET needs 5mC at CpGs for assignment')
    wanted = [x.strip() for x in a.mods.split(',') if x.strip()]
    for mod in wanted:
        if mod in ('5hmCG', '5hmCH') and not has_h:
            errors.append('%s requested but the BAM has no 5hmC (C+h) calls '
                          '(basecalled with a 5mCG or 5mC model)' % mod)
        if mod in ('5mCH', '5hmCH') and not all_context:
            errors.append('%s requested but the BAM only has CpG-context calls '
                          '(basecalled with a 5mCG or 5mCG_5hmCG model)' % mod)
    for e in errors:
        sys.stderr.write('ERROR: %s\n' % e)
    print(kind if has_m else 'unsupported')
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
