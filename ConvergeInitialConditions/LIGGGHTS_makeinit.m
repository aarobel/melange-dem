%===========================================================
%====== EXAMPLE 1: input file for the CONVERGENCE example
%===========================================================
% log-normal FSD, no bonds
N = 1600;
rmean = 4.5;
rsigma = 0.8;

params = struct('N'   ,N,   ...  % number of floes (-)
                'arat',2.0,     ...  % model domain's aspect ratio Lx/Ly (-)
                'A',0.3,        ...  % ice concentration (-)
                'rho',910.0,    ...  % ice density (kg/m3)
                'h',100.0,        ...  % ice thickness (m)
                'Vwrite',0);         % write the 'Velocities' section?
rdist = struct('type','logn','rmean',rmean,'sigma',rsigma,'rmax',500,'rmin',30);
Vdist = struct('type','norm','Vrndmeanx',0.0,'Vrndstdx',1e-10, ...
                             'Vrndmeany',0.0,'Vrndstdy',1e-10); 
parbnd = struct('ibnd'   ,0,    ...  % whether bond-related info should be written to the header
                                ...  %  (0=no bond info; 1=only space for bonds left, but no. of bonds = 0; 
                                ...  %   2=bonds created between all grains separated by less than dmax) 
                'nextra',5,     ...  % No. of "extra bond per atom"
                'dmax',5.0,     ...  % 
                'nbondtypes',2, ...  % No. of bond types
                'bondtyperatio',...
                     [0.8 0.2], ...  % proportions of various types of bonds
                'bondstoremove',...
                     0.0        ...  % proportion of bonds to be removed from the total
                );
[L,r,Rx,Ry,Vx,Vy] = LIGGGHTSinit(params,rdist,Vdist,parbnd,['conv_N',int2str(N),'_mean',num2str(rmean),'_std',num2str(rsigma),'_Ly7.init']);
