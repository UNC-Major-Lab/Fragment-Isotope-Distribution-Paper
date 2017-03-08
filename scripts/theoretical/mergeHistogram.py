#!/usr/bin/env python
import sys
import os
import re
from collections import defaultdict
from math import floor
from math import isnan

root_dir = sys.argv[1]
prefix = sys.argv[2]
do_iso = sys.argv[3]

if do_iso == 'F':
    comp2bin2count = defaultdict(dict)

    for f in os.listdir(root_dir):
        fp = root_dir+"/"+f
        if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
            infile = open(fp)
            for line in infile:
                [comp, bin, count] = line.strip().split("\t")
                count = int(count)
                if not comp2bin2count[comp].has_key(bin):
                    comp2bin2count[comp][bin] = 0
                comp2bin2count[comp][bin]+=count
            close(fp)

    for comp in comp2bin2count:
        for bin in comp2bin2count[comp]:
            print "\t".join([comp, str(bin), str(comp2bin2count[comp][bin])])

else:
    comp2iso2bin2count = defaultdict(dict)

    for f in os.listdir(root_dir):
        fp = root_dir+"/"+f
        if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):

            infile = open(fp)
            for line in infile:
                [comp, bin, iso, count] = line.strip().split("\t")
                count = int(count)
                if not comp2iso2bin2count[comp].has_key(iso):
                    comp2iso2bin2count[comp][iso] = defaultdict(int)
                comp2iso2bin2count[comp][iso][bin]+=count
            close(fp)

    for comp in comp2iso2bin2count:
        for iso in comp2iso2bin2count[comp]:
            for bin in comp2iso2bin2count[comp][iso]:
                print "\t".join([comp, str(bin), iso, str(comp2iso2bin2count[comp][iso][bin])])