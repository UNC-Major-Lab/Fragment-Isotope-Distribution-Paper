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


method2stats = defaultdict(list)
fileCount = 0;

for f in os.listdir(root_dir):
    fp = root_dir+"/"+f
    if os.path.isfile(fp) and ".out" in f and f.startswith(prefix):
        infile = open(fp)
        fileCount+=1

        for line in infile:
            if do_iso == "F":
                [count, mean, min_v, q1, median, q3, max_v, method, stat] = line.strip().split("\t")
                label = method +"\t"+stat
            else:
                [count, mean, min_v, q1, median, q3, max_v, iso, method, stat] = line.strip().split("\t")
            label = iso + "\t" + method +"\t"+stat

            count = int(count)
            mean = float(mean)
            min_v = float(min_v)
            q1 = float(q1)
            median = float(median)
            q3 = float(q3)
            max_v = float(max_v)


            if not method2stats.has_key(label):
                method2stats[label] = [count, count * mean, min_v, q1, median, q3, max_v]
            else:
                method2stats[label][0] += count
                method2stats[label][1] += count * mean
                method2stats[label][2] = min(min_v,  method2stats[label][2])
                method2stats[label][3] += q1
                method2stats[label][4] += median
                method2stats[label][5] += q3
                method2stats[label][6] = max(max_v,  method2stats[label][6])

        infile.close()

for method in method2stats:
    mean = method2stats[method][1] / method2stats[method][0]
    min_v =  method2stats[method][2]
    q1 =  method2stats[method][3] / fileCount
    median =  method2stats[method][4] / fileCount
    q3 =  method2stats[method][5] / fileCount
    max_v =  method2stats[method][6]

    print mean + "\t" + min_v + "\t" + q1 + "\t" + median + "\t" + q3 + "\t" + max_v
