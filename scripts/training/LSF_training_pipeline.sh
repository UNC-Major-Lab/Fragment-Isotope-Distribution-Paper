#!/bin/sh

if ($1 == 1) then
    bsub < LSF_get_sulfur_distribution.sh
    bsub < LSF_create_average_training_data.sh -w 'ended('LSF_get_sulfur_distribution.sh')'
    bsub < LSF_create_training_data.sh
else if ($1 == 2) then
    bsub < LSF_combine_models.sh
endif