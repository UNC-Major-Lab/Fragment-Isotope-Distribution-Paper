#!/bin/csh

bsub < LSF_compare_to_theoretical.sh > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out
bsub < LSF_merge_theoretical.sh -w 'numended('$JOBID',>=300)'