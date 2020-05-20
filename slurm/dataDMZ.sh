#!/bin/bash

# Script to create scratch directory and download data to it from gcloud
# This should be implemented in the Science DMZ nodes
# e.g. dtn01-dmz.chpc.utah.edu

# Define start time (for clocking execution speed)
t_start=`date +%s`

# Purge old modules and load required ones
module purge
module load rclone

# Define location of input data directories
DATADIR="gcloud:CHPC/IceBridge-raw/wais-central/20111109/"

# Define location of input density file
RHOFILE="gcloud:CHPC/flight-density/rho_20111109.csv"

# Define and create scratch directory
SCRDIR="/scratch/general/lustre/u1046484/"
mkdir -p $SCRDIR

# Transfer input data to scratch
rclone copyto $RHOFILE $SCRDIR/rho_data.csv
rclone copy $DATADIR $SCRDIR/Data/20111109/ --transfers=16 --drive-chunk-size=32768

# Define end time and calculate execution time
t_end=`date +%s`
runtime=$((t_end-t_start))
echo "Total execution time:"
echo $runtime
