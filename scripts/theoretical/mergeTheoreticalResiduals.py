#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
from math import floor
from math import isnan

root_dir = sys.argv[1]
bin_size = float(sys.argv[2])
prefix = sys.argv[3]

bin2count = defaultdict(int)

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and prefix in f:
        infile = open(fp)
        for line in infile:
            residual = line.strip()
            bin = round(float(residual)/bin_size)*bin_size
            bin2count[bin]+=1

for bin in bin2count:
    print "\t".join([str(bin), str(bin2count[bin])])
