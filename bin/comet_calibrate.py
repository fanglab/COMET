#!/usr/bin/env python3
"""Pass-number-aware ML cutoff calibration for PacBio HiFi (Methods 2.4).

Input : ML histogram from comet_pacbio_scan.py (np, counts of ML 0..255)
Output: pass-number table, one row per np:
        np, tau0, meth(tau0), unmeth(tau0), tau_np, meth(tau_np), unmeth(tau_np)
        tau_np = 256 means calls from that pass number are excluded.

With --no-adjust every pass number keeps tau0 (fixed cutoff).
"""
import sys
import math
import argparse

import numpy as np


def count(h, tau):
    return int(h[tau:].sum()), int(h[:256 - tau].sum())  # ML >= tau, ML <= 255 - tau


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('hist')
    p.add_argument('out')
    p.add_argument('--cutoff', type=int, default=180)
    p.add_argument('--no-adjust', action='store_true')
    a = p.parse_args()
    cutoff = a.cutoff

    fd = {}
    for line in open(a.hist):
        f = line.split('\t')
        fd[int(f[0])] = np.array([int(x) for x in f[1:]], dtype=np.int64)

    adjust = not a.no_adjust
    if adjust:
        # sample-level reference methylation from high-accuracy reads (15 <= np <= 20)
        m = sum(count(v, cutoff)[0] for k, v in fd.items() if 15 <= k <= 20)
        u = sum(count(v, cutoff)[1] for k, v in fd.items() if 15 <= k <= 20)
        if m + u == 0:
            sys.exit('ERROR: no confident calls from reads with 15-20 passes; '
                     'cannot calibrate. Re-run with --no-pn-calibration.')
        glb = m / (m + u)

    fd.pop(1, None)
    fd.pop(2, None)
    with open(a.out, 'w') as out:
        for n in sorted(fd):
            ct = fd[n]
            if ct.sum() == 0:
                continue
            meth, unmeth = count(ct, cutoff)
            total = meth + unmeth
            row = [n, cutoff, meth, unmeth]
            if not adjust:
                row.extend([cutoff, meth, unmeth])
            elif total == 0:
                row.extend([256, 0, 0])
            else:
                ml_c = meth / total
                if n >= 15 or glb - ml_c <= 0.005:
                    row.extend([cutoff, meth, unmeth])
                else:
                    # raise tau one step at a time; keep the candidate window
                    # |M_n(tau) - g| within 0.5% of g, stop right after it
                    flag = 0
                    ml2_d = {}
                    for n2 in np.arange(1, 255 - cutoff):
                        cutoff2 = cutoff + int(n2)
                        meth2, unmeth2 = count(ct, cutoff2)
                        total2 = meth2 + unmeth2
                        if total2 > 0 and math.isclose(meth2 / total2, glb, rel_tol=0.005):
                            flag = 1
                            ml2_d[meth2 / total2] = [meth2, unmeth2, cutoff2]
                        elif flag == 1:
                            break
                    if ml2_d:
                        best = min(ml2_d, key=lambda k: abs(k - glb))
                        meth2, unmeth2, cutoff2 = ml2_d[best]
                        row.extend([cutoff2, meth2, unmeth2])
                    else:
                        row.extend([256, 0, 0])
            out.write('\t'.join(str(i) for i in row) + '\n')


if __name__ == '__main__':
    main()
