#!/bin/bash
#SBATCH --job-name=sicalve
#SBATCH --output=liggghts_calve.out
#SBATCH --error=liggghts_calve.err
#SBATCH --time=1:00:00
#SBATCH --partition=sandyb
#SBATCH --account=pi-abbot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32000

module load liggghts/3.1+intelmpi-5.0+intel-15.0

mpirun -np 4 lmp_intelmpi -echo screen < in.springcalve
