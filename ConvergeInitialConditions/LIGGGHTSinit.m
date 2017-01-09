function [L,r,Rx,Ry,Vx,Vy] = ...
                         LIGGGHTSinit(params,rdist,Vdist,parbnd,outfile)
% LIGGGHTSinit - writes an input file for LAMMPS/LIGGGHTS 
%                with floe radii and initial floe positions and velocities.
%
% Author: Agnieszka Herman, IOUG (agnieszka.herman@ug.edu.pl)
%
N   = params.N;
rng('shuffle')
%========================================================================
%=== Calculate N and floe radii:
%========================================================================
if isequal(rdist.type,'constant')
    r = rdist.rmean*ones(N,1);
elseif isequal(rdist.type,'unif')
    r = random('unif',rdist.rmin,rdist.rmax,N,1);
elseif isequal(rdist.type,'logn')
    r = random('logn',rdist.rmean,rdist.sigma,N,1);
    for(q=1:N);while(r(q)>rdist.rmax | r(q)<rdist.rmin);r(q) = random('logn',rdist.rmean,rdist.sigma,1,1);end;end
elseif isequal(rdist.type,'gp')
    r = random('gp',rdist.k,rdist.sigma,rdist.rmean,N,1);
elseif isequal(rdist.type,'bidisp')
    n1 = round((rdist.rmean-rdist.r2)*N/(rdist.r1-rdist.r2));
    r = [rdist.r1*ones(n1,1); rdist.r2*ones(N-n1,1)];
% elseif isequal(rdist.type,'glv')
%     Smax = A*L^2/pi;
%     r = [];
%     x = rdist.xmin:0.1:rdist.xmax;
%     P = x.^(-1-rdist.alpha).*exp((1-rdist.alpha)./x);    
%     Ntmp = 1e3;
%     while sum(r.^2)<Smax
%         r = rdist.rmean*randpdf(P,x,[Ntmp 1]);
%         Ntmp = Ntmp + 1e3;
%     end
%     N = find(cumsum(r.^2)-Smax>0,1,'first');
%     r = r(1:N);
elseif isequal(rdist.type,'powerlaw')
    r = ((rdist.alpha*N+1)./(rdist.alpha*(1:N)'+1)).^(1/rdist.alpha);
    r = r.*exp(-r.^2*rdist.alpha/1e5); % correction so that the largest
                                        % floes are not too large
    r = r/mean(r)*rdist.rmean;
end
r = sort(r,'descend');
Ly = sqrt(pi*sum(r.^2)/params.A/params.arat);
Lx = params.arat*Ly;
L = [Lx Ly];
%========================================================================
%=== Initialize floes' velocities:
%========================================================================
if params.Vwrite > 0
    if isequal(Vdist.type,'norm')
        Vx = random('norm',Vdist.Vrndmeanx,Vdist.Vrndstdx,N,1);
        Vy = random('norm',Vdist.Vrndmeany,Vdist.Vrndstdy,N,1);
    end
else
    Vx = [];
    Vy = [];
end
%========================================================================
%=== Initialize floe positions:
%========================================================================
Rx = random('unif',-Lx/2,Lx/2,N,1);
Ry = random('unif',-Ly/2,Ly/2,N,1);
for ii = 2:N
    disp(ii);
    overlap = 1;
    while overlap > 0
        overlap = 0;
        for j = 1:ii-1
            RXIJ = Rx(ii) - Rx(j);
            RYIJ = Ry(ii) - Ry(j);
            if (RXIJ - Lx*round(RXIJ/Lx))^2 + ...
               (RYIJ - Ly*round(RYIJ/Ly))^2 < (r(ii)+r(j))^2
                overlap = 1;
            end
        end
        if overlap==1
            Rx(ii) = random('unif',-Lx/2,Lx/2,1,1);
            Ry(ii) = random('unif',-Ly/2,Ly/2,1,1);
        end
    end   
end
%========================================================================
%=== Initialize bonds:
%========================================================================
if parbnd.ibnd == 2
    bondlist = [];
    nextra = 0;
    for ii = 1:N-1
        ipartner = ii+1:N;
        RXIJ = Rx(ii) - Rx(ipartner);
        RYIJ = Ry(ii) - Ry(ipartner);
        % bonds are not created around a periodic boundary:
        dist = sqrt(RXIJ.^2 + RYIJ.^2) - (r(ii)+r(ipartner));
        ind = find(dist <= parbnd.dmax);
        nextra = max(nextra,length(ind));
        ibnd = [ii*ones(length(ind),1) ipartner(ind)'];
        bondlist = [bondlist; ibnd];
    end
    if parbnd.bondstoremove > 0
        if parbnd.bondstoremove > 1
            error('parbnd.bondstoremove must be lower than 1');
        else
            ind = randperm(size(bondlist,1),round(parbnd.bondstoremove*size(bondlist,1)));
            bondlist(ind,:) = [];
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
    end
end
%========================================================================
%=== Write the LIGGGHTS input file(s):
%========================================================================
if parbnd.ibnd == 0
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',1,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',-Lx/2,Lx/2,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-0.5,0.5,' zlo zhi');
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
        fprintf(fid,'%6.0f%2.0f%9.3f%9.3f%9.3f%15.3f%15.3f%15.3f\n', ...
                             n,1,2*r(n),params.h,params.rho,Rx(n),Ry(n),0);
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f%9.3f%9.3f%9.3f%4.1f%4.1f%4.1f\n',n,Vx(n),Vy(n),0,0,0,0);
        end
    end
    fclose(fid);
elseif parbnd.ibnd == 1
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',1,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',-Lx/2,Lx/2,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-0.5,0.5,' zlo zhi');
    fprintf(fid,'%5.0f%s\n\n',0,' bonds');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nbondtypes,' bond types');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nextra,' extra bond per atom');    
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
        fprintf(fid,'%6.0f%2.0f%15.3f%15.3f%15.3f%12.3f%9.3f%9.3f%2.0f\n', ...
                             n,1,Rx(n),Ry(n),0,2*r(n),params.h,params.rho,1);
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f%9.3f%9.3f%9.3f%4.1f%4.1f%4.1f\n',n,Vx(n),Vy(n),0,0,0,0);
        end
    end
    fclose(fid);
elseif parbnd.ibnd == 2
    fid = fopen(outfile,'w');
    fprintf(fid,'%s\n\n',' Input file for LIGGGHTS, created with Matlab function LIGGGHTSinit.m');
    fprintf(fid,'%5.0f%s\n\n',N,' atoms');
    fprintf(fid,'%5.0f%s\n\n',1,' atom types');
    fprintf(fid,'%15.3f%15.3f%s\n',-Lx/2,Lx/2,' xlo xhi');
    fprintf(fid,'%15.3f%15.3f%s\n',-Ly/2,Ly/2,' ylo yhi');
    fprintf(fid,'%15.3f%15.3f%s\n\n',-0.5,0.5,' zlo zhi');
    fprintf(fid,'%5.0f%s\n\n',size(bondlist,1),' bonds');    
    fprintf(fid,'%5.0f%s\n\n',parbnd.nbondtypes,' bond types');    
    fprintf(fid,'%5.0f%s\n\n',max(parbnd.nextra,nextra),' extra bond per atom');    
    fprintf(fid,'%s\n\n',' Atoms');
    for n = 1:N
        fprintf(fid,'%6.0f%2.0f%15.3f%15.3f%15.3f%12.3f%9.3f%9.3f%2.0f\n', ...
                             n,1,Rx(n),Ry(n),0,2*r(n),params.h,params.rho,1);
    end
    if params.Vwrite > 0
        fprintf(fid,'\n%s\n\n',' Velocities');
        for n = 1:N
            fprintf(fid,'%6.0f%9.3f%9.3f%9.3f%4.1f%4.1f%4.1f\n',n,Vx(n),Vy(n),0,0,0,0);
        end
    end
    fprintf(fid,'\n%s\n\n',' Bonds');
    for n = 1:size(bondlist,1)
        fprintf(fid,'%8.0f%3.0f%8.0f%8.0f\n',n,bondtype(n),bondlist(n,1),bondlist(n,2));
    end
    fclose(fid);    
end
