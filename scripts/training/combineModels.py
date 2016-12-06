#!/usr/bin/env python
import sys
import os
import re

root_dir = sys.argv[1]

print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
print "<models>\n"

for d in os.listdir(root_dir):
    if not os.path.isfile(d) and os.path.exists(root_dir+d+"/spline/model/"):
        for f in os.listdir(root_dir+d+"/spline/model/"):
            fp = root_dir+d+"/spline/model/"+f
            if os.path.isfile(fp):
                infile = open(fp)
                for line in infile:
                    print line.rstrip()

print "</models>"
