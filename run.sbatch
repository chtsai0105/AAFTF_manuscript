#!/bin/bash

#SBATCH --job-name=aaftf_benchmark
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/aaftf_benchmark.%j.log

pixi run nextflow run workflow/main.nf -params-file config/params_test.yaml -profile pixi,slurm -resume
