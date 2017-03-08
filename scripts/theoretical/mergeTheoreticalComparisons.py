#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
from math import floor
from math import isnan

root_dir = sys.argv[1]
num_bins = float(sys.argv[2])
prefix = sys.argv[3]
job_id = int(sys.argv[4])
num_jobs = int(sys.argv[5])

comp2iso2bin2count = dict()

max_job = num_jobs * job_id
min_job = max_job - num_jobs

max_val = -float("inf")
min_val = float("inf")

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        job = int(f.split("_")[-1].split(".")[0])
        if job <= min_job or job > max_job: continue

        infile = open(fp)
        for line in infile:
            [score, iso, comp] = line.strip().split("\t")
            score = float(score)
            max_val = max(max_val, score)
            min_val = min(min_val, score)

bin_size = (max_val-min_val) / num_bins

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        job = int(f.split("_")[-1].split(".")[0])
        if job <= min_job or job > max_job: continue

        infile = open(fp)
        for line in infile:
            [score, iso, comp] = line.strip().split("\t")
            if not comp2iso2bin2count.has_key(comp):
                comp2iso2bin2count[comp] = dict()
            if not comp2iso2bin2count[comp].has_key(iso):
                comp2iso2bin2count[comp][iso] = defaultdict(int)
            bin = round(float(score)/bin_size)*bin_size
            comp2iso2bin2count[comp][iso][bin]+=1

for comp in comp2iso2bin2count:
    for iso in comp2iso2bin2count[comp]:
        for bin in comp2iso2bin2count[comp][iso]:
            print "\t".join([comp, str(bin), iso, str(comp2iso2bin2count[comp][iso][bin])])
