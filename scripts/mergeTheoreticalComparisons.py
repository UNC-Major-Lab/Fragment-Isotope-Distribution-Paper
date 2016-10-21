#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
from math import floor
from math import isnan

root_dir = sys.argv[1]

comp2iso2bin2count = dict()

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f:
        infile = open(fp)
        for line in infile:
            [pcc, man, chi, pep_mass, frag_mass, iso, num_s, num_cs, comp] = line.strip().split("\t")
            if not comp2iso2bin2count.has_key(comp):
                comp2iso2bin2count[comp] = dict()
            if not comp2iso2bin2count[comp].has_key(iso):
                comp2iso2bin2count[comp][iso] = defaultdict(int)
            bin = round(float(man),2)
            comp2iso2bin2count[comp][iso][bin]+=1

for comp in comp2iso2bin2count:
    for iso in comp2iso2bin2count[comp]:
        for bin in comp2iso2bin2count[comp][iso]:
            print "\t".join([comp, str(bin), iso, str(comp2iso2bin2count[comp][iso][bin])])
