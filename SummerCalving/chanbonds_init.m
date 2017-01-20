res_orig = readdump(['./res/compress_gaussnarrow500_lgberg.all.7920000']);

atomdata_orig = assembleresults_fix(res_orig);
params.rho = 910; params.h = 100; params.Vwrite = 0;

btr=1;
for Lx = 12e3:1e3:16e3;

narrowing_length_scale = 5e2;
narrowing_x_center = 0;
narrowing_height = 1e3;

wall_function = @(x) narrowing_height*exp(-((x-narrowing_x_center)/narrowing_length_scale).^2);

%for btr = [1]
      for pf = 0.700:0.005:.770 
        r = atomdata_orig.r;
        Rx = atomdata_orig.x;
        Ry = atomdata_orig.y;
        N  = length(r);
        mlg_idx = find(abs(Ry)<2.5e3-wall_function(Rx) & Rx>min(Rx)+4e2 & Rx<min(Rx)+10.4e3);  
        packingfrac_orig = pi.*sum(r(mlg_idx).^2)./((max((r(mlg_idx)+Rx(mlg_idx)))-min((-r(mlg_idx)+Rx(mlg_idx))))*(max((r(mlg_idx)+Ry(mlg_idx)))-min((-r(mlg_idx)+Ry(mlg_idx)))));

        r_orig = r;

        for shrinkfactor = linspace(1,0.8,1e4)
              r(mlg_idx) = shrinkfactor.*r_orig(mlg_idx);
              packingfrac = pi.*sum(r(mlg_idx).^2)./((max((r(mlg_idx)+Rx(mlg_idx)))-min((-r(mlg_idx)+Rx(mlg_idx))))*(max((r(mlg_idx)+Ry(mlg_idx)))-min((-r(mlg_idx)+Ry(mlg_idx))))); 
              packingfrac*100
              if(packingfrac<pf);break;end
        end
	idx = find(atomdata_orig.x>min(atomdata_orig.x)+Lx & abs(atomdata_orig.y)<2.5e3);	
	idxr = idx(randperm(length(idx),floor(0.4*length(idx))));        


	r(idxr) = [];
	Rx(idxr) = [];
	Ry(idxr) = [];

	atomdata = atomdata_orig;
        atomdata.r = [r;10;10];
        atomdata.x = [Rx;max(Rx)-500;max(Rx)-500];
        atomdata.y = [Ry;-1e3;1e3];
        %atomdata.r = [r];
        %atomdata.x = [Rx];
        %atomdata.y = [Ry];
	
	atomdata.x = atomdata.x-(min(atomdata.x)+range(atomdata.x)/2); %reset center of domain to zero
	atomdata.y = atomdata.y-(min(atomdata.y)+range(atomdata.y)/2); %reset mean of center to zero

	atomdata.L = [max(atomdata.x+atomdata.r)-min(atomdata.x-atomdata.r),max(atomdata.y+atomdata.r)-min(atomdata.y-atomdata.r)];
	
	for dmax = [10]

		parbnd = struct('ibnd'   ,3,    ...  % whether bond-related info should be written to the header
                                ...  %  (0=no bond info; 1=only space for bonds left, but no. of bonds = 0; 
                                ...  %   2=bonds created between all grains separated by less than dmax) 
                'nextra',100,     ...  % No. of "extra bond per atom"
                'dmax',dmax,     ...  % 
                'nbondtypes',1, ...  % No. of bond types
                'bondtyperatio',...
                     [1.0 0.0], ...  % proportions of various types of bonds
                'bondstoremove',...
                     btr,        ...  % proportion of bonds to be removed from the total
                'term_width',...
                     210,        ...  %terminus width
                'wall_width',...
                     10,        ...  %sidewall width
                'wall_radius',...
		     atomdata.r(find(atomdata.y==min(atomdata.y),1)));
        
	[L,r,Rx,Ry,Vx,Vy] = LIGGGHTSinit_fromdump(atomdata,params,parbnd,['./init/fjordbond_lgberg_gn1000_btr' int2str(btr*100) '_dmax' int2str(dmax) '_pf' int2str(pf*1e3) '_Lx' int2str(Lx/1e3) '.init']);
	
	end
   end
end
