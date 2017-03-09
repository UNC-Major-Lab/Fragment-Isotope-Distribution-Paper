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
job_id = int(sys.argv[4])-1
num_jobs = int(sys.argv[5])

comp2bin2count = defaultdict(dict)

comp2max_val = dict()
comp2min_val = dict()
comp2bin_size = defaultdict(float)



for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and prefix in f:
        job = int(f.split("_")[-1].split(".")[0])
        if job % num_jobs != job_id: continue

        infile = open(fp)
        for line in infile:
            [residual, comp] = line.strip().split("\t")
            residual = float(residual)

            if not comp2max_val.has_key(comp):
                comp2max_val[comp] = -float("inf")
                comp2min_val[comp] = float("inf")

            comp2max_val[comp] = max(comp2max_val[comp], residual)
            comp2min_val[comp] = min(comp2min_val[comp], residual)
        infile.close()

for comp in comp2max_val:
    comp2bin_size[comp] = (comp2max_val[comp]-comp2min_val[comp]) / num_bins

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and prefix in f:
        job = int(f.split("_")[-1].split(".")[0])
        if job % num_jobs != job_id: continue

        infile = open(fp)
        for line in infile:
            [residual, comp] = line.strip().split("\t")
            bin = round(float(residual)/comp2bin_size[comp])*comp2bin_size[comp]
            if not comp2bin2count[comp].has_key(bin):
                comp2bin2count[comp][bin]=0
            comp2bin2count[comp][bin]+=1
        infile.close()

for comp in comp2bin2count:
    for bin in comp2bin2count[comp]:
        print "\t".join([comp, str(bin), str(comp2bin2count[comp][bin])])
