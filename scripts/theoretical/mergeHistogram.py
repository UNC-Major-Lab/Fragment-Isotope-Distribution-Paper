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

comp2bin2count = defaultdict(dict)

max_job = num_jobs * job_id
min_job = max_job - num_jobs

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        job = int(f.split("_")[-1].split(".")[0])
        if job <= min_job or job > max_job: continue

        infile = open(fp)
        for line in infile:
            [comp, bin, count] = line.strip().split("\t")
            count = int(count)
            if not comp2bin2count[comp].has_key(bin):
                comp2bin2count[comp][bin] = 0
            comp2bin2count[comp][bin]+=count

for comp in comp2bin2count:
    for bin in comp2bin2count[comp]:
        print "\t".join([comp, str(bin), str(comp2bin2count[comp][bin])])
