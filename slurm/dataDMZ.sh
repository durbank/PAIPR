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
DATADIR="gcloud:CHPC/IceBridge-raw/ALL_flightlines/"

# Define filtering file to download subset of echograms
FILTER="filter-file.txt"

# Define location of input density file
RHODIR="gcloud:CHPC/flight-density/"

echo "Creating scratch directory"

# Define and create scratch directory
# SCRDIR="/scratch/general/lustre/u1046484/"
SCRDIR="/scratch/general/nfs1/u1046484/"
mkdir -p $SCRDIR

# Transfer input density data to scratch
echo "Transfering density data"
rclone copy $RHODIR $SCRDIR/rho_data/ --progress

# Transfer input radar data (based on filtering if it exists)
if [ -e $FILTER ]
then
  echo "Transfering filtered echogram data"
  rclone copy $DATADIR $SCRDIR/Data/ --filter-from=$FILTER --transfers=16 --drive-chunk-size=32768 --progress
else
  echo "Transfering all echogram data"
  rclone copy $DATADIR $SCRDIR/Data/ --transfers=16 --drive-chunk-size=32768 --progress
fi

# Define end time and calculate execution time
t_end=`date +%s`
runtime=$((t_end-t_start))
echo "Total execution time:"
echo $runtime
