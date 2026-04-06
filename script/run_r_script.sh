#!/bin/bash

#SBATCH --job-name=sc_rna_analysis
#SBATCH --account=def-cottenie
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=250GB
#SBATCH --error=%x_%j_error.txt
#SBATCH --output=%x_%j_output.txt

module load r/4.5.0

#allow user to input their R script
INPUT=${1}

#line to run inputted R script
Rscript $INPUT



