#!/bin/bash
set -e  # stop on first error

cd 'Population Analyses'

Rscript --no-save --slave '00_Source_All.R'

cd ..

cd 'Simulation Study'

bash run_scripts.sh