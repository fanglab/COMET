#!/usr/bin/env python3
"""Scan a PacBio HiFi modBAM once (one worker per autosome).

Outputs
  <prefix>.reads.tsv    read_id, chrom, ref_start, np   (reads to be scored)
  <prefix>.ml_hist.tsv  np, then counts of ML values 0..255 (for calibration)

Read selection follows the paper pipeline:
  calibration histogram : primary, mapped (-F 0x904), autosomal, length >= min_len
  scored reads          : same, plus flag in {0, 16} and np > 2
"""
import re
import sys
import argparse
from multiprocessing import Pool

import numpy as np
import pysam

AUTOSOME = re.compile(r'^chr([1-9]|1[0-9]|2[0-2])$')


def scan_chrom(args):
    bam_file, chrom, min_len = args
    hist = {}
    reads = []
    no_np = 0
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam.fetch(chrom):
            if read.flag & 0x904:
                continue
            if read.query_length < min_len:
                continue
            if not read.has_tag('np'):
                no_np += 1
                continue
            npc = int(read.get_tag('np'))
            if read.has_tag('ML'):
                ml = np.frombuffer(read.get_tag('ML'), dtype=np.uint8)
                counts = np.bincount(ml, minlength=256)
            else:
                counts = np.zeros(256, dtype=np.int64)
            if npc in hist:
                hist[npc] += counts
            else:
                hist[npc] = counts.astype(np.int64)
            if (read.flag == 0 or read.flag == 16) and npc > 2:
                reads.append((read.query_name, chrom, read.reference_start, npc))
    return chrom, hist, reads, no_np


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bam')
    p.add_argument('out_prefix')
    p.add_argument('--min-len', type=int, default=7500)
    p.add_argument('--threads', type=int, default=1)
    a = p.parse_args()

    with pysam.AlignmentFile(a.bam, 'rb') as bam:
        chroms = [c for c in bam.references if AUTOSOME.match(c)]
    if not chroms:
        sys.exit('ERROR: no autosomes named chr1..chr22 in BAM header')

    jobs = [(a.bam, c, a.min_len) for c in chroms]
    with Pool(max(1, min(a.threads, len(jobs)))) as pool:
        results = pool.map(scan_chrom, jobs)

    hist = {}
    no_np = 0
    seen = set()
    n_reads = 0
    with open(a.out_prefix + '.reads.tsv', 'w') as f:
        for chrom, h, reads, nn in results:  # pool.map keeps chromosome order
            no_np += nn
            for k, v in h.items():
                hist[k] = hist[k] + v if k in hist else v
            for rid, c, pos, npc in reads:
                if rid in seen:
                    continue
                seen.add(rid)
                n_reads += 1
                f.write('%s\t%s\t%d\t%d\n' % (rid, c, pos, npc))

    with open(a.out_prefix + '.ml_hist.tsv', 'w') as f:
        for k in sorted(hist):
            f.write(str(k) + '\t' + '\t'.join(str(int(x)) for x in hist[k]) + '\n')

    if no_np:
        sys.stderr.write('WARNING: %d reads without an np tag were skipped\n' % no_np)
    sys.stderr.write('Scored-read candidates: %d\n' % n_reads)
    if n_reads == 0:
        sys.exit('ERROR: no reads passed filters (length >= %d, np > 2, primary, autosomal)' % a.min_len)


if __name__ == '__main__':
    main()
