#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --mem=48G

#SBATCH -J GB_spinodal_Ia
#SBATCH -o GB_spinodal_Ia.out
#SBATCH -e GB_spinodal_Ia.out

mkdir output
./spinodal_gb
