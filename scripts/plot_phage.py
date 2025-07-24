#!/usr/bin/env python3
# coding: utf-8

import bacchus.plot as bcp
import cooler
import sys

cool_file = sys.argv[1]
res = int(sys.argv[2])
out_file = sys.argv[3]

mat = cooler.Cooler(f'{cool_file}').matrix(balance=True)[:]
chrom_starts = bcp.get_chrom_start(cool_file, res)
bcp.contact_map(mat, binning=res, vmax=95, chrom_starts=chrom_starts, out_file=out_file)
