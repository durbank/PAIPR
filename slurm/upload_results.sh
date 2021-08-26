#!/bin/bash

# Template script for uploading PAIPR-generated results to gcloud

# Assign remote directory to which to upload
CLOUD_DIR="gcloud:CHPC/PAIPR-results"

# Assign path to source directory
SRC_DIR=$PWD

# Get name of enclosing direcotry
dir_name="$(basename $SRC_DIR)"

# Load the rclone module
module load rclone

# Define start time (for clocking execution speed)
t_start=`date +%s`

# Transfer sbatch source files and results to gcloud
echo "Transferring files to gcloud..."
rclone copy $SRC_DIR $CLOUD_DIR/$dir_name --transfers=8 --progress
echo "Tranfer complete"

# Define end time and calculate execution time
t_end=`date +%s`
runtime=$((t_end-t_start))
echo "Total execution time:"
echo $runtime

