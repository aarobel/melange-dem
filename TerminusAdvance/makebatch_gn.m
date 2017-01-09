function makebatch(nt,iter,run,hsi,str,lfc)

runname = [run '_h' int2str(10*hsi) '_str' num2str(str/1e6) 'e6_lfc' int2str(1000*lfc)];
shrt = ['_gn250_h' int2str(10*hsi) '_str' num2str(str/1e6) 'e6_lfc' int2str(1000*lfc)];

fid = fopen(['batch_' runname '.sh'],'w');

fprintf(fid,'%s\n','#!/bin/bash');
fprintf(fid,'%s\n',['#SBATCH --job-name=bcgnh' int2str(10*hsi) 's' num2str(str/1e6) 'l' int2str(lfc*1000)]);
fprintf(fid,'%s\n',['#SBATCH --output=liggghts' shrt '.out']);
fprintf(fid,'%s\n',['#SBATCH --error=liggghts' shrt '.err']);
fprintf(fid,'%s\n','#SBATCH --time=36:00:00');
fprintf(fid,'%s\n','#SBATCH --partition=sandyb');
fprintf(fid,'%s\n','#SBATCH --account=pi-abbot');
fprintf(fid,'%s\n','#SBATCH --nodes=1');
fprintf(fid,'%s\n','#SBATCH --ntasks-per-node=16');
fprintf(fid,'%s\n\n','#SBATCH --mem=32000');

fprintf(fid,'%s\n','module load liggghts/3.1+intelmpi-5.0+intel-15.0');
fprintf(fid,'%s\n','module load matlab');

fprintf(fid,'%s\n',['matlab -nodesktop -nosplash -r "chanbonds_fromstop_gn(' sprintf('\''') run sprintf('\''') ',$1,$2,' num2str(hsi) ',' int2str(str) ',' num2str(lfc) ');makeLAMMPSrun_gn($2,' sprintf('\''') run sprintf('\''') ',' num2str(hsi) ',' int2str(str) ',' num2str(lfc) ');quit;"']);
fprintf(fid,'%s\n',['mpirun -np 16 lmp_intelmpi -echo screen < in.' runname]);

fclose(fid);
