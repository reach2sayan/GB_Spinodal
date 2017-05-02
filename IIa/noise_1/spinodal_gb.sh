#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=48G

#SBATCH -J GB_spinodal_IIa
#SBATCH -o GB_spinodal_IIa.out
#SBATCH -e GB_spinodal_IIa.out

mkdir output
./spinodal_gb