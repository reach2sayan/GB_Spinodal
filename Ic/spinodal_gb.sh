#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --mem=48G

#SBATCH -J GB_spinodal_Ic
#SBATCH -o GB_spinodal_Ic.out
#SBATCH -e GB_spinodal_Ic.out

mkdir output
./spinodal_gb