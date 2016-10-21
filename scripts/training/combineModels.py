#!/usr/bin/env python
import sys
import os
import re

schema_infile = open(sys.argv[1])
root_dir = sys.argv[2]

for line in schema_infile:
    print line.rstrip()

for d in os.listdir(root_dir):
    if not os.path.isfile(d) and os.path.exists(root_dir+d+"/spline/model/"):
        for f in os.listdir(root_dir+d+"/spline/model/"):
            fp = root_dir+d+"/spline/model/"+f
            if os.path.isfile(fp):
                infile = open(fp)
                for line in infile:
                    print line.rstrip()
print "</models>\n</class>"
