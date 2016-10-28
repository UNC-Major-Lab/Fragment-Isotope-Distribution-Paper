bsub < compare_to_theoretical.sh > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out
bsub < merge_theoretical.sh -w 'numended('$JOBID',>=300)'