function [L,r,Rx,Ry,Vx,Vy] = ...
                         LIGGGHTSinit_boncreate(atomdata,bondlist,params,parbnd,outfile)
% LIGGGHTSinit_fromdump - writes an input file for LAMMPS/LIGGGHTS 
%                         based on data read from a dump file
%
% Author: Agnieszka Herman, IOUG (agnieszka.herman@ug.edu.pl)
%
N  = length(atomdata.r);
r = atomdata.r;
Lx = atomdata.L(1);
Ly = atomdata.L(2);
xlo = -Lx/2;
xhi = Lx/2;

L = [xhi-xlo,Ly];

Rx = atomdata.x;
Ry = atomdata.y;


    rxmin = min(Rx)+parbnd.term_width;
    rymin = -2.5e3;
    rymax = 2.5e3;

    newbondlist = [];
    nextra = 0;
    rlist = [];
   if(parbnd.nolf==1)
        for ii = 1:N-1
                %if(Rx(ii)>rxmin && Ry(ii)>rymin+wall_function(Rx(ii)) && Ry(ii)<rymax-wall_function(Rx(ii)))
           if((r(ii)~=parbnd.term_radius & r(ii)~=parbnd.wall_radius))
                ipartner = ii+1:N;
                RXIJ = Rx(ii) - Rx(ipartner);
                RYIJ = Ry(ii) - Ry(ipartner);
                % bonds are not created around a periodic boundary:
                dist = sqrt(RXIJ.^2 + RYIJ.^2) - (r(ii)+r(ipartner));
                ind = find(dist <= parbnd.dmax & (r(ipartner)~=parbnd.term_radius & r(ipartner)~=parbnd.wall_radius));
                nextra = max(nextra,length(ind));
                ibnd = [ii*ones(length(ind),1) ipartner(ind)'];
                newbondlist = [newbondlist; ibnd];
                rlist = [rlist; r(ii)*ones(length(ind),1) r(ipartner(ind)')];
           end
        end
   else
        for ii = 1:N-1
                %if(Rx(ii)>rxmin && Ry(ii)>rymin+wall_function(Rx(ii)) && Ry(ii)<rymax-wall_function(Rx(ii)))
           if(r(ii)~=parbnd.term_radius)
                ipartner = ii+1:N;
                RXIJ = Rx(ii) - Rx(ipartner);
                RYIJ = Ry(ii) - Ry(ipartner);
                % bonds are not created around a periodic boundary:
                dist = sqrt(RXIJ.^2 + RYIJ.^2) - (r(ii)+r(ipartner));
                ind = find(dist <= parbnd.dmax & r(ipartner)~=parbnd.term_radius);
                nextra = max(nextra,length(ind));
                ibnd = [ii*ones(length(ind),1) ipartner(ind)'];
                newbondlist = [newbondlist; ibnd];
                rlist = [rlist; r(ii)*ones(length(ind),1) r(ipartner(ind)')];
           end
        end
   end
for jj = 1:length(newbondlist)
   if(~((newbondlist(jj,1)==bondlist(:,1) & newbondlist(jj,2)==bondlist(:,2)) | ...
        (newbondlist(jj,2)==bondlist(:,1) & newbondlist(jj,1)==bondlist(:,2))))
	bondlist = [bondlist;newbondlist(jj,:)];
   end
end

    if parbnd.nbondtypes>1
        Nb = size(bondlist,1);
        bondtype = zeros(Nb,1);
        nbtot = 0;
        for ii = 1:parbnd.nbondtypes-1
            nbt = round(Nb*parbnd.bondtyperatio(ii));
            bondtype(nbtot+1:nbtot+nbt) = ii;
            nbtot = nbtot + nbt;
        end
        bondtype(nbtot+1:end) = parbnd.nbondtypes;
        bondtype = bondtype(randperm(Nb));
    else
        Nb = size(bondlist,1);
        bondtype = ones(Nb,1);
    end


if params.Vwrite > 0
    Vx = atomdata.vx;
    Vy = atomdata.vy;
    omegax = atomdata.omegax;
    omegay = atomdata.omegay;
    omegaz = atomdata.omegaz;
else
    Vx = [];
    Vy = [];
end

%========================================================================
%=== Write the LIGGGHTS input file(s):
%========================================================================
if parbnd.ibnd == 0
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit_fromdump.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',2,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',xlo,xhi,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-100,100,' zlo zhi');
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
     if(r(n)~=parbnd.wall_radius)
          fprintf(fid,'%6.0f%2.0f%12.3f%9.3f%9.3f%15.3f%15.3f%15.3f\n', ...
                             n,1,2*r(n),2*r(n),params.rho,Rx(n),Ry(n),0); 
        else
           fprintf(fid,'%6.0f%2.0f%12.3f%9.3f%9.3f%15.3f%15.3f%15.3f\n', ...
                             n,2,2*r(n),2*r(n),params.rho,Rx(n),Ry(n),0); 
        end  
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f%9.3f%9.3f%9.3f%12.5f%12.5f%12.5f\n',n,Vx(n),Vy(n),0,omegax(n),omegay(n),omegaz(n));
        end
    end
    fclose(fid);
elseif parbnd.ibnd == 1
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit_fromdump.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',2,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',xlo,xhi,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-100,100,' zlo zhi');
    fprintf(fid,'%5.0f%s\n\n',0,' bonds');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nbondtypes,' bond types');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nextra,' extra bond per atom');    
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
        if(r(n)==parbnd.wall_radius)
           fprintf(fid,'%6.0f%2.0f%15.3f%15.3f%15.3f%12.3f%9.3f%9.3f%2.0f\n', ...
                             n,2,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
        else if(r(n)==parbnd.term_radius)
           fprintf(fid,'%6.0f%2.0f%15.3f%15.3f%15.3f%12.3f%9.3f%9.3f%2.0f\n', ...
                             n,3,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
        else
           fprintf(fid,'%6.0f%2.0f%15.3f%15.3f%15.3f%12.3f%9.3f%9.3f%2.0f\n', ...
                             n,1,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
        end;end     
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f%9.3f%9.3f%9.3f%12.5f%12.5f%12.5f\n',n,Vx(n),Vy(n),0,omegax(n),omegay(n),omegaz(n));
        end
    end
    fclose(fid);
elseif parbnd.ibnd == 2 || parbnd.ibnd == 3
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit_fromdump.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',2,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',xlo,xhi,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-100,100,' zlo zhi');
    fprintf(fid,'%5.0f%s\n\n',size(bondlist,1),' bonds');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nbondtypes,' bond types');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nextra,' extra bond per atom');    
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
	if(r(n)==parbnd.wall_radius)
           fprintf(fid,'%6.0f%2.0f%15.6f%15.6f%15.6f%12.3f%12.3f%9.3f%2.0f\n', ...
                             n,2,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
	else if(r(n)==parbnd.term_radius)
	   fprintf(fid,'%6.0f%2.0f%15.6f%15.6f%15.6f%12.3f%12.3f%9.3f%2.0f\n', ...
                             n,1,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
	else
           fprintf(fid,'%6.0f%2.0f%15.6f%15.6f%15.6f%12.3f%12.3f%9.3f%2.0f\n', ...
                             n,1,Rx(n),Ry(n),0,2*r(n),2*r(n),params.rho,1);
	end;end
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f %9.8f %9.8f %9.8f %12.8f %12.8f %12.8f\n',n,Vx(n),Vy(n),0,omegax(n),omegay(n),omegaz(n));
        end
    end
    fprintf(fid,'\n%s\n\n',' Bonds');
    for n = 1:size(bondlist,1)
        fprintf(fid,'%8.0f%3.0f%8.0f%8.0f\n',n,bondtype(n),bondlist(n,1),bondlist(n,2));
    end
    fclose(fid);    
end
