#!/bin/sh

bsub < LSF_speed_test.sh

bsub < LSF_compare_to_theoretical.sh > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out

bsub < LSF_merge_theoretical2.sh -w 'numended('$JOBID',>=100)'