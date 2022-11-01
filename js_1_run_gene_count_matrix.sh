#!/bin/bash
#SBATCH --mem=80960 --cpus-per-task=20
#SBATCH --time=48:00:00

source activate r_deseq

R CMD BATCH run_gene_count_matrix.R
