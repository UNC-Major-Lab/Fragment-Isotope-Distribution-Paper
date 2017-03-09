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

comp2iso2bin2count = dict()

comp2max_val = dict()
comp2min_val = dict()
comp2bin_size = defaultdict(float)

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        job = int(f.split("_")[-1].split(".")[0])
        if job % num_jobs != job_id: continue

        infile = open(fp)
        for line in infile:
            [score, iso, comp] = line.strip().split("\t")
            score = float(score)


            if not comp2max_val.has_key(comp):
                comp2max_val[comp] = -float("inf")
                comp2min_val[comp] = float("inf")

            comp2max_val[comp] = max(comp2max_val[comp], score)
            comp2min_val[comp] = min(comp2min_val[comp], score)

        infile.close()

for comp in comp2max_val:
    comp2bin_size[comp] = (comp2max_val[comp]-comp2min_val[comp]) / num_bins

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        job = int(f.split("_")[-1].split(".")[0])
        if job % num_jobs != job_id: continue

        infile = open(fp)
        for line in infile:
            [score, iso, comp] = line.strip().split("\t")
            if not comp2iso2bin2count.has_key(comp):
                comp2iso2bin2count[comp] = dict()
            if not comp2iso2bin2count[comp].has_key(iso):
                comp2iso2bin2count[comp][iso] = defaultdict(int)
            bin = round(float(score)/comp2bin_size[comp])*comp2bin_size[comp]
            comp2iso2bin2count[comp][iso][bin]+=1
        infile.close()

for comp in comp2iso2bin2count:
    for iso in comp2iso2bin2count[comp]:
        for bin in comp2iso2bin2count[comp][iso]:
            print "\t".join([comp, str(bin), iso, str(comp2iso2bin2count[comp][iso][bin])])
