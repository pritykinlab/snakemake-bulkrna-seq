#!/bin/bash
### create a tmp dir for the pairtools
source activate snakemake_bulkRNA
### run make snake in the background
nohup bash snakemake.sh &
