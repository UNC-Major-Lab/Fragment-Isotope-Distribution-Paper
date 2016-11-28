#!/bin/csh

bsub < LSF_get_sulfur_distribution.sh > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`

bsub < LSF_create_training_data.sh > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out
bsub < LSF_combine_models.sh -w 'numended('$JOBID',>=16)'

bsub < LSF_create_average_training_data.sh -w 'ended('LSF_get_sulfur_distribution.sh')' > out
JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out