function makeLAMMPSrun(iter,run,hsi,str,lfc)

runname = [run '_h' int2str(10*hsi) '_str' num2str(str/1e6) 'e6_lfc' int2str(lfc*1000)]

fid = fopen(['in.' runname],'w');
fprintf(fid,'%s\n\n','#====== Important Variables');
fprintf(fid,'%s\n',['variable        runID   string "' runname '"']);
fprintf(fid,'%s\n','variable        dt      equal  0.002');
fprintf(fid,'%s\n','variable        dtsnap  equal  60*60/${dt}');
fprintf(fid,'%s\n','variable        dtthrm  equal  10*60/${dt}');
fprintf(fid,'%s\n','variable        trun    equal  5*24*3600/${dt}');
fprintf(fid,'%s\n','variable	termvel equal  40/(24*3600)');
fprintf(fid,'%s\n','variable        bondupd equal  24*60*60/${dt}');
fprintf(fid,'%s\n',['variable	iter	equal  ',int2str(iter)]);
fprintf(fid,'%s\n\n','variable        bergradius         equal  1e3');


fprintf(fid,'%s\n\n','#===== Initialization:');
fprintf(fid,'%s\n','units           si');
fprintf(fid,'%s\n','dimension	2');
fprintf(fid,'%s\n','newton		off off');
fprintf(fid,'%s\n','boundary	p p p');
fprintf(fid,'%s\n\n','processors	8 2 1');

fprintf(fid,'%s\n\n','#===== Atom definitions:');
fprintf(fid,'%s\n','atom_style	hybrid disk bond/gran/disk');
fprintf(fid,'%s\n','atom_modify	map array');
fprintf(fid,'%s\n','pair_style      gran model hertz/stiffness/disk tangential_damping on');
fprintf(fid,'%s\n','bond_style	gran/disk');
fprintf(fid,'%s\n\n','special_bonds	lj/coul 0 1 1 extra 40');

if(iter==1)
        fprintf(fid,'%s\n',['read_data       ./init/fjordbond_gn250_lgberg_nlf_mbtr0_dmax10.init']);
else
        fprintf(fid,'%s\n',['read_data       ./res/' runname '.restart']);
end

fprintf(fid,'%s\n','pair_coeff	* *');
fprintf(fid,'%s\n','fix             m1 all property/global youngsModulus peratomtype 9.e9 9.9e3');
fprintf(fid,'%s\n','fix             m2 all property/global poissonsRatio peratomtype 0.33 0.33');
fprintf(fid,'%s\n','fix             m3 all property/global coefficientRestitution peratomtypepair 2 0.9 0.9 0.9 0.9');
fprintf(fid,'%s\n','fix             m4 all property/global coefficientFriction peratomtypepair 2 0.7 0.7 0.7 0.7');
fprintf(fid,'%s\n','fix             m5 all property/global kn peratomtypepair 2 6.7e8 6.7e8 6.7e8 6.7e8');
fprintf(fid,'%s\n','fix             m6 all property/global kt peratomtypepair 2 2.7e8 2.7e8 2.7e8 2.7e8');
fprintf(fid,'%s\n','fix             m7 all property/global gamman peratomtypepair 2 1.0 1.0 1.0 1.0');
fprintf(fid,'%s\n\n','fix             m8 all property/global gammat peratomtypepair 2 0.5 0.5 0.5 0.5');
%fprintf(fid,'%s\n\n','fix           m9 all viscous 2e5');

fprintf(fid,'%s\n\n','#===== Bond properties:');
fprintf(fid,'%s\n','variable        E     equal 9.e9 # Young (elastic) modulus');
fprintf(fid,'%s\n','variable        kn2ks equal 1.5  # normal to shear stiffness ratio');
fprintf(fid,'%s\n',['variable        sigmacmax equal ' int2str(str)]);
fprintf(fid,'%s\n',['variable        sigmatmax equal ' int2str(str)]);
fprintf(fid,'%s\n',['variable        taumax    equal ' int2str(str)]);
fprintf(fid,'%s\n\n',['bond_coeff      1 1.0 ' num2str(hsi) ' 0.0 0.01 ${E} ${kn2ks} 1 ${sigmacmax} ${sigmatmax} ${taumax}']);

fprintf(fid,'%s\n\n','#===== General settings:');
fprintf(fid,'%s\n','neighbor	1.0 nsq');
fprintf(fid,'%s\n','neigh_modify	delay 0 page 10000 one 300');
fprintf(fid,'%s\n','timestep	${dt}');
fprintf(fid,'%s\n\n','communicate	single cutoff 6000.0 vel yes');

fprintf(fid,'%s\n\n','#===== Forcing-related settings:');
fprintf(fid,'%s\n','fix             1  all nve/disk');
fprintf(fid,'%s\n','fix             f2b all property/global rhoWater scalar 1025.0');
fprintf(fid,'%s\n','fix		ff2 all property/global rhoAir	 scalar 1.27');
fprintf(fid,'%s\n\n','fix             f3  all enforce2d');

fprintf(fid,'%s\n\n','#======= define glacier terminus subregion:');

fprintf(fid,'%s\n','variable        yleft        equal bound(all,ymin)');
fprintf(fid,'%s\n','variable        ycleft       equal ${yleft}+210');
fprintf(fid,'%s\n','variable        yright       equal bound(all,ymax)');
fprintf(fid,'%s\n','variable        ycright      equal ${yright}-210');
fprintf(fid,'%s\n','variable        xbot         equal -10100');
fprintf(fid,'%s\n','variable        xtbot        equal -10010+(${termvel}*${trun}*${dt}*(${iter}-1))');
fprintf(fid,'%s\n','variable        xtop         equal bound(all,xmax)');
fprintf(fid,'%s\n','variable        xtermbot     equal ${xtbot}+220');
fprintf(fid,'%s\n','variable        xtermtop     equal ${xtop}-210');
%fprintf(fid,'%s\n','variable        xlengthpos   equal ${xtermtop}-${xtermbot}');
%fprintf(fid,'%s\n','variable        xlengthneg   equal ${xtermbot}-${xtermtop}');
%fprintf(fid,'%s\n','variable        ywalldiffneg equal -1*${ywalldiffpos}');
%fprintf(fid,'%s\n','variable        ycleftmin    equal ${ycleft}+${ywalldiffpos}');
%fprintf(fid,'%s\n\n','variable        ycrightmin   equal ${ycright}-${ywalldiffpos}');

fprintf(fid,'%s\n\n','group           chanwalls type 2');

fprintf(fid,'%s\n','region          termb block ${xtbot} ${xtermbot} ${ycleft} ${ycright} -100 100 units box');
fprintf(fid,'%s\n\n','group           term region termb');

fprintf(fid,'%s\n','region          botl block ${xbot} ${xtermbot} ${yleft} ${ycleft} -100 100 units box');
fprintf(fid,'%s\n','region          botr block ${xbot} ${xtermbot} ${ycright} ${yright} -100 100 units box');
fprintf(fid,'%s\n','region          topl block ${xtermtop} ${xtop} ${yleft} ${ycleft} -100 100 units box');
fprintf(fid,'%s\n','region          topr block ${xtermtop} ${xtop} ${ycright} ${yright} -100 100 units box');
fprintf(fid,'%s\n','region          extrabits union 4 botl botr topl topr');
fprintf(fid,'%s\n','group           ebt region extrabits');
fprintf(fid,'%s\n','group           cfrz union chanwalls ebt');

fprintf(fid,'%s\n','group           notmelange union cfrz term');
fprintf(fid,'%s\n','group           melange subtract all notmelange');
fprintf(fid,'%s\n','group           noterm  subtract all term');

fprintf(fid,'%s\n\n','#===== Output:');
fprintf(fid,'%s\n','compute         2b all stress/atom pair');
fprintf(fid,'%s\n','compute         2a all stress/atom bond');
fprintf(fid,'%s\n','compute         bstr all reduce sum c_2a[1] c_2a[2] c_2a[4]');
fprintf(fid,'%s\n','compute         pstr all reduce sum c_2b[1] c_2b[2] c_2b[4]');
fprintf(fid,'%s\n','compute         kerot all erotate/disk');
fprintf(fid,'%s\n','compute         4a all pair/gran/local id force');
fprintf(fid,'%s\n','compute         4b all property/local batom1 batom2');
fprintf(fid,'%s\n','compute         4d all bond/gran/disk/local length thickness forcen forcet torquetz');
fprintf(fid,'%s\n','compute         5a all pressure thermo_temp pair');
fprintf(fid,'%s\n','compute         5b all pressure thermo_temp bond');
fprintf(fid,'%s\n\n','fix             fm all ave/atom 1 ${dtsnap} ${dtsnap} vx vy fx fy c_2b[1] c_2b[2] c_2b[4] c_2a[1] c_2a[2] c_2a[4]');
fprintf(fid,'%s\n','thermo          ${dtthrm}');
fprintf(fid,'%s\n','thermo_style    custom time c_5a[1] c_5a[2] c_5a[4] c_5b[1] c_5b[2] c_5b[4] c_bstr[1] c_bstr[2] c_bstr[3] c_pstr[1] c_pstr[2] c_pstr[3]');
fprintf(fid,'%s\n','thermo_style    custom time c_kerot ke pxx pyy pxy lx ly');
fprintf(fid,'%s\n',['dump            1 melange custom 20000000 ./res/${runID}.iter${iter}.melange.const id radius mass']);
fprintf(fid,'%s\n',['dump            2 melange custom ${dtsnap} ./res/${runID}.iter${iter}.melange.* id radius xs ys vx vy omegaz fx fy tqz c_2b[1] c_2b[2] c_2b[4]']);

fprintf(fid,'%s\n',['dump            4 term custom 20000000 ./res/${runID}.iter${iter}.terminus.const id radius mass']);
fprintf(fid,'%s\n',['dump            5 term custom ${dtsnap} ./res/${runID}.iter${iter}.terminus.* id radius xs ys vx vy ']);
fprintf(fid,'%s\n',['dump            6 all local  ${dtsnap} ./res/${runID}.iter${iter}.all.*.pairs c_4a[1] c_4a[2] c_4a[3] c_4a[4] c_4a[5]']);
fprintf(fid,'%s\n',['dump            7 all local  ${dtsnap} ./res/${runID}.iter${iter}.all.*.bonds c_4b[1] c_4b[2] c_4d[1] c_4d[2] c_4d[3] c_4d[4] c_4d[5]']);
fprintf(fid,'%s\n',['dump            8 notmelange custom  ${dtsnap} ./res/${runID}.iter${iter}.notmelange.* id radius xs ys vx vy omegaz fx fy tqz c_2b[1] c_2b[2] c_2b[4]']);
fprintf(fid,'%s\n\n',['dump          9 all custom ${dtsnap} ./res/${runID}.iter${iter}.all.* id radius xs ys vx vy omegaz fx fy tqz c_2b[1] c_2b[2] c_2b[4]']);


fprintf(fid,'%s\n\n','#===== Run:');
fprintf(fid,'%s\n','#Start with no calving at terminus');
fprintf(fid,'%s\n','fix             frz    cfrz  freeze');
fprintf(fid,'%s\n\n','fix             ftermove term move linear ${termvel} 0 0 units box');

fprintf(fid,'%s\n','run             ${trun}');

fclose(fid);
