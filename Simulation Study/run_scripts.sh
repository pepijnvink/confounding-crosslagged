#!/bin/bash
set -e  # stop on first error
mkdir -p output
mkdir -p figures

for i in main bothdir highconf
do
  mkdir -p output/$i
  mkdir -p output/$i/simdat
  mkdir -p output/$i/analyzed_dat
  mkdir -p output/$i/results
  mkdir -p figures/$i
  mkdir -p figures/$i/png
  mkdir -p figures/$i/png/x
  mkdir -p figures/$i/png/y
  mkdir -p figures/$i/power
  mkdir -p figures/$i/svg
  mkdir -p figures/$i/svg/x
  mkdir -p figures/$i/svg/y
  for j in 250 500 1000 10000
  do
      mkdir -p output/$i/analyzed_dat/N$j
      mkdir -p output/$i/results/N$j
      mkdir -p output/$i/simdat/N$j
  done
done

for i in main bothdir highconf
do
  for j in 250 500 1000 10000
  do
		Rscript --no-save --slave 'scripts/$i/00_$i_multiplication_tables.R'
		Rscript --no-save --slave 'scripts/$i/01_$i_simulate.R' $j 2500
		Rscript --no-save --slave 'scripts/$i/02_$i_analyses.R' $j
		Rscript --no-save --slave 'scripts/$i/03_$i_plot_results.R' $j x
		Rscript --no-save --slave 'scripts/$i/03_$i_plot_results.R' $j y
  done
  Rscript --no-save --slave 'scripts/$i/03_$i_plot_power.R' x
  Rscript --no-save --slave 'scripts/$i/03_$i_plot_power.R' y
done
