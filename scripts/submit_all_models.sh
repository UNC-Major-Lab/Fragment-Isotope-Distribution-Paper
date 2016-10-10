#!/bin/csh 

set max_mass = $1

./submit_models.sh 0 0 0 0 $max_mass
./submit_models.sh 1 0 0 0 $max_mass
./submit_models.sh 1 1 0 0 $max_mass
./submit_models.sh 2 0 0 0 $max_mass
./submit_models.sh 2 1 0 0 $max_mass
./submit_models.sh 2 2 0 0 $max_mass
./submit_models.sh 3 0 0 0 $max_mass
./submit_models.sh 3 1 0 0 $max_mass
./submit_models.sh 3 2 0 0 $max_mass
./submit_models.sh 3 3 0 0 $max_mass
./submit_models.sh 4 0 0 0 $max_mass
./submit_models.sh 4 1 0 0 $max_mass
./submit_models.sh 4 2 0 0 $max_mass
./submit_models.sh 5 0 0 0 $max_mass
./submit_models.sh 5 1 0 0 $max_mass
./submit_models.sh 6 0 0 0 $max_mass
