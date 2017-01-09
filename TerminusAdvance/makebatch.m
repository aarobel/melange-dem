function makebatch(nt,iter,run,hsi,str,lfc)

runname = [run '_h' int2str(10*hsi) '_str' num2str(str/1e5) 'e5_lfc' int2str(1000*lfc)];
shrt = ['_h' int2str(10*hsi) '_str' num2str(str/1e5) 'e5_lfc' int2str(1000*lfc)];

fid = fopen(['batch_' runname '.sh'],'w');

fprintf(fid,'%s\n','#!/bin/bash');
fprintf(fid,'%s\n',['#SBATCH --job-name=bch' int2str(10*hsi) 'sd' num2str(str/1e5)]);
fprintf(fid,'%s\n',['#SBATCH --output=liggghts' shrt '.out']);
fprintf(fid,'%s\n',['#SBATCH --error=liggghts' shrt '.err']);
fprintf(fid,'%s\n','#SBATCH --time=24:00:00');
fprintf(fid,'%s\n','#SBATCH --partition=sandyb');
fprintf(fid,'%s\n','#SBATCH --account=pi-abbot');
fprintf(fid,'%s\n','#SBATCH --nodes=1');
fprintf(fid,'%s\n','#SBATCH --ntasks-per-node=8');
fprintf(fid,'%s\n\n','#SBATCH --mem=32000');

fprintf(fid,'%s\n','module load liggghts/3.1+intelmpi-5.0+intel-15.0');
fprintf(fid,'%s\n','module load matlab');

fprintf(fid,'%s\n',['matlab -nodesktop -nosplash -r "chanbonds_fromstop(' sprintf('\''') run sprintf('\''') ',$1,$2,' num2str(hsi) ',' int2str(str) ',' num2str(lfc) ');makeLAMMPSrun($2,' sprintf('\''') run sprintf('\''') ',' num2str(hsi) ',' int2str(str) ',' num2str(lfc) ');quit;"']);
fprintf(fid,'%s\n',['mpirun -np 8 lmp_intelmpi -echo screen < in.' runname]);

fclose(fid);
