#!/bin/bash
#SBATCH --job-name=convinit
#SBATCH --output=liggghts_init.out
#SBATCH --error=liggghts_init.err
#SBATCH --time=4:00:00
#SBATCH --partition=viz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load matlab
matlab -nodisplay < LIGGGHTS_makeinit.m 

module load liggghts/3.1+intelmpi-5.0+intel-15.0
mpirun -np 1 lmp_intelmpi -echo screen < in.conv
