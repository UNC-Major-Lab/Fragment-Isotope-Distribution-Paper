#!/bin/sh

source ../config.sh

bsub < LSF_get_sulfur_distribution.sh

bsub < LSF_create_training_data.sh > out

JOBID=`head -1 out | sed 's/.*<\\([0-9]*\\)>.*/\\1/'`
rm out

TRAINING_JOBS = 7

bsub < LSF_combine_models.sh  -w 'numended('$JOBID',>='$TRAINING_JOBS')'
