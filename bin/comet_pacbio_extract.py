#!/usr/bin/env python3
"""Per-read CpG 5mC calls from a PacBio HiFi modBAM for one chunk of reads (Methods 2.2).

Input : chunk read list (read_id, chrom, ref_start, np) and pass-number table
Output: chrom, pos, pos+1, ML, np, strand, read_id   (one line per CpG call)

C+m calls are taken at the C; G-m calls are shifted one base to the C of the CpG.
A call is kept only if the read has 'CG' there, the base is an aligned match, and
the call is confident under the read's pass-number cutoff (ML >= tau or ML <= 255 - tau).
Only the genomic span of the chunk is fetched from the BAM.
"""
import argparse
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bam')
    p.add_argument('chunk_reads')
    p.add_argument('pn_table')
    p.add_argument('out')
    a = p.parse_args()

    span = {}
    targets = set()
    order = []
    for line in open(a.chunk_reads):
        rid, chrom, pos = line.split('\t')[:3]
        pos = int(pos)
        targets.add(rid)
        if chrom not in span:
            span[chrom] = [pos, pos]
            order.append(chrom)
        else:
            span[chrom][0] = min(span[chrom][0], pos)
            span[chrom][1] = max(span[chrom][1], pos)

    pn_df = pd.read_csv(a.pn_table, sep='\t', header=None, index_col=0)[[4]]
    pn_d = pn_df.to_dict()[4]

    bam = pysam.AlignmentFile(a.bam, 'rb')
    with open(a.out, 'w') as f:
        for chrom in order:
            start, end = span[chrom]
            for read in bam.fetch(chrom, start, end + 1):
                if read.reference_start < start or read.query_name not in targets:
                    continue
                if not (read.flag == 0 or read.flag == 16):
                    continue
                npc = read.get_tag('np')
                if npc <= 2:
                    continue
                targets.discard(read.query_name)  # one record per read
                seq = read.query_sequence
                mod_bases = read.modified_bases or {}
                # 5mC is tagged 'C+m' (qpos at the C of CG) or 'G-m' (qpos at the G, C at qpos-1)
                calls = []
                for (canon_base, mstrand, mod_code), pairs in mod_bases.items():
                    if canon_base != 'C' or mod_code != 'm':
                        continue
                    offset = 0 if mstrand == 0 else -1
                    calls.extend((qpos + offset, ml) for qpos, ml in pairs)
                if not calls:
                    continue
                match = np.array(read.get_aligned_pairs(matches_only=True))
                if match.size == 0:
                    continue
                match_d = dict(zip(match[:, 0], match[:, 1].astype(int)))
                strand = '-' if read.is_reverse else '+'
                tau = pn_d.get(npc, 256)
                for cpos, ml in calls:
                    if cpos not in match_d or seq[cpos:cpos + 2] != 'CG':
                        continue
                    if ml >= tau or ml <= 255 - tau:
                        ref_pos = match_d[cpos]
                        f.write('%s\t%d\t%d\t%d\t%d\t%s\t%s\n' % (
                            read.reference_name, ref_pos, ref_pos + 1, ml, npc,
                            strand, read.query_name))
    bam.close()


if __name__ == '__main__':
    main()
