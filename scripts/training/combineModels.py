#!/usr/bin/env python
import sys
import os
import re
from datetime import datetime

def iterateAndPrint(root_dir):
    for d in os.listdir(root_dir):
        if not os.path.isfile(d) and os.path.exists(root_dir+d+"/spline/model/"):
            for f in os.listdir(root_dir+d+"/spline/model/"):
                fp = root_dir+d+"/spline/model/"+f
                if os.path.isfile(fp):
                    infile = open(fp)
                    for line in infile:
                        print line.rstrip()

root_dir = sys.argv[1]
max_depth = sys.argv[2]
max_sulfur = sys.argv[3]

now = datetime.utcnow();

print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
print "<models maxIsotope=\"" + max_depth + "\" maxSulfur=\"" + max_sulfur + "\" createdDate=\"" + now.isoformat(' ') + "\">"

iterateAndPrint(root_dir);

print "</models>"
