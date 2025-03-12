#!/bin/bash
# This script installs and sets up the Souporcell pipeline

# Load Singularity module (check available versions with 'module avail singularity')
module load singularity/3.7.1

# Create a new directory for Souporcell
mkdir souporcell
cd souporcell || exit 1

# Download Singularity image for Souporcell
singularity pull shub://wheaton5/souporcell

# Download Souporcell pipeline
wget https://raw.githubusercontent.com/wheaton5/souporcell/master/souporcell_pipeline.py

# Download Demuxafy Singularity image and MD5 checksum
wget -O Demuxafy.sif 'https://www.dropbox.com/scl/fi/g0cuyjwomdavom6u6kb2v/Demuxafy.sif?rlkey=xfey1agg371jo4lubsljfavkh&'
wget -O Demuxafy.sif.md5 'https://www.dropbox.com/scl/fi/bk3p2k2440un6sb6psijn/Demuxafy.sif.md5?rlkey=x3vl8ejpfhjsrvmjanwzkxty9'

# Download Assign_Indiv_by_Geno.R script
wget https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs/blob/main/scripts/Assign_Indiv_by_Geno.R

echo "Installation and setup completed successfully."
