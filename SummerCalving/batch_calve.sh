#!/bin/bash
#SBATCH --job-name=bndcalve
#SBATCH --output=liggghts_cv.out
#SBATCH --error=liggghts_cv.err
#SBATCH --time=2:00:00
#SBATCH --partition=viz
#SBATCH --account=pi-abbot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

module load liggghts/3.1+intelmpi-5.0+intel-15.0

mpirun -np 4 lmp_intelmpi -echo screen < in.summercalve
