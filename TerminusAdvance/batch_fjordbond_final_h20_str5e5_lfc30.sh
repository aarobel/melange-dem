#!/bin/bash
#SBATCH --job-name=bch20sd5
#SBATCH --output=liggghts_h20_str5e5_lfc30.out
#SBATCH --error=liggghts_h20_str5e5_lfc30.err
#SBATCH --time=24:00:00
#SBATCH --partition=sandyb
#SBATCH --account=pi-abbot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=32000

module load liggghts/3.1+intelmpi-5.0+intel-15.0
module load matlab
matlab -nodesktop -nosplash -r "chanbonds_fromstop('fjordbond_final',$1,$2,2,500000,0.03);makeLAMMPSrun($2,'fjordbond_final',2,500000,0.03);quit;"
mpirun -np 8 lmp_intelmpi -echo screen < in.fjordbond_final_h20_str5e5_lfc30
