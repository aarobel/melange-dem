function chanbonds_fromstop_nlfgn(runfilename,nt,iter,hsi,str,lfc)
if(iter~=1)
runname = [runfilename '_h' int2str(10*hsi) '_str' num2str(str/1e6) 'e6_lfc' int2str(lfc*1000)]
runfile = [runname '.iter' int2str(iter-1) '.all.' int2str(nt)]

res_orig = readdump(['./res/' runfile ]);
bondump = readdump(['./res/' runfile '.bonds']);
bondlist = bondump.entries{1}(:,1:2);

atomdata_orig = assembleresults_fix(res_orig);
params.rho = 910; params.h = 100; params.Vwrite = 1;

atomdata = atomdata_orig;
%atomdata.L = [max(atomdata.x+atomdata.r)-min(atomdata.x-atomdata.r),max(atomdata.y+atomdata.r)-min(atomdata.y-atomdata.r)];
r_mid = atomdata.r(abs(atomdata.y)<1e3);
x_mid = atomdata.x(abs(atomdata.y)<1e3);
dmax = 10;

Lx = 16;

narrowing_length_scale = 5e2;
narrowing_x_center = 0;
narrowing_height = 0;

wall_function = @(x) narrowing_height*exp(-((x-narrowing_x_center)/narrowing_length_scale).^2);

parbnd = struct('ibnd'   ,3,    ...  % whether bond-related info should be written to the header
                                ...  %  (0=no bond info; 1=only space for bonds left, but no. of bonds = 0; 
                                ...  %   2=bonds created between all grains separated by less than dmax) 
                'nextra',40,     ...  % No. of "extra bond per atom"
                'dmax',10,     ...  % 
                'nbondtypes',1, ...  % No. of bond types
                'bondtyperatio',...
                     [1.0 0.0], ...  % proportions of various types of bonds
                'bondstoremove',...
                     0,        ...  % proportion of bonds to be removed from the total
                'term_width',...
                     210,        ...  %terminus width
                'wall_width',...
                     10,        ...  %sidewall width
                'wall_radius',...
		     atomdata.r(find(atomdata.y==min(atomdata.y),1)),...
		'term_radius',...
		     atomdata.r(1490),...
                'narrow_height',narrowing_height,...
                'nolf',1);
		[L,r,Rx,Ry,Vx,Vy] = LIGGGHTSinit_bondcreate(atomdata,bondlist,params,parbnd,['./res/' runname '.restart']);

end
end
