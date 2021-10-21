#!/bin/bash

# Development script for testing bash script functionality locally

# Define current working directory
startDIR=$(pwd)

SRCFILE="/home/durbank/Documents/Research/Antarctica/PAIPR/src/scripts/process_local.m"
# directory with the Matlab functions required to run src script
SRCDIR="/home/durbank/Documents/Research/Antarctica/PAIPR/src/"
# Path to density results file
RHOFILE="/media/durbank/WARP/Research/Antarctica/Data/rho-stats/Data/cent_interior/Data/densities/rho_20111109.csv"
# directory with the data input files
DATADIR="/media/durbank/WARP/Research/Antarctica/Data/IceBridge/WAIS-central/20111109/"
# directory to output results
OUTDIR=$startDIR/tmp_results/
# OUTDIR="/home/durbank/Documents/MATLAB/scratch/tmp_results/"
# Assign scratch directory to process data
SCRDIR="/home/durbank/Documents/scratch/"

# create results directory in home space
mkdir -p $OUTDIR

# make scratch dir, copy input data there and cd to it
mkdir -p $SCRDIR
cp $SRCFILE $SCRDIR
cp -r $SRCDIR $SCRDIR/src
cp -r $DATADIR $SCRDIR/Data
cp $RHOFILE $SCRDIR/rho_data.csv
cd $SCRDIR

# Create seperate directory within scratch to store outputs
mkdir Outputs

# run matlab program with variable options
matlab -nodisplay -r "process_local" -logfile Outputs/matlab_batch.log

# copy results from scratch to the result dir
cp -r Outputs/ $OUTDIR

# Return to initial working directory
cd $startDIR

# Clear the scratch directory
rm -rf $SCRDIR
