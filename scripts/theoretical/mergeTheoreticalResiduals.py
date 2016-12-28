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
job_id = int(sys.argv[4])
num_jobs = int(sys.argv[5])

bin2count = defaultdict(int)

max_job = num_jobs * job_id
min_job = max_job - num_jobs

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and prefix in f:
        job = f.split("_")[-1].split(".")[0]
        if job <= min_job or job > max_job: continue

        infile = open(fp)
        for line in infile:
            residual = line.strip()
            bin = round(float(residual)/bin_size)*bin_size
            bin2count[bin]+=1

for bin in bin2count:
    print "\t".join([str(bin), str(bin2count[bin])])
