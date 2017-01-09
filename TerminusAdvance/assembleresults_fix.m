function atomdata = assembleresults_fix(var)
% assembleresults -- assembles atom data into a single structure  
% that may be used, e.g., as input to LIGGGHTSinit_fromdump or 
% to plotfloepositions.
%
% Author: Agnieszka Herman, IOUG (agnieszka.herman@ug.edu.pl)
%
% k = findposition(const,'id');
% [id,ind] = sort(squeeze(const.data(:,k)),'ascend');
% k = findposition(const,'radius');
% r  = squeeze(const.data(:,k));
% r = r(ind);
% L(1) = var.xb(2)-var.xb(1);
% L(2) = var.yb(2)-var.yb(1);
k = findposition(var,'id');
[id,ind] = sort(var.data(:,k),'ascend');
k = findposition(var,'radius');
r  = var.data(ind,k);
% r = r(ind);
L(1) = var.xb(2)-var.xb(1);
L(2) = var.yb(2)-var.yb(1);
k = findposition(var,'xs');
x = var.data(ind,k);
k = findposition(var,'ys');
y = var.data(ind,k);
x = x*L(1) - L(1)/2;
y = y*L(2) - L(2)/2;
atomdata = struct('id',id,'r',r,'L',L,'x',x,'y',y);
atomdata.vx = readposition(var,'vx',ind);
atomdata.vy = readposition(var,'vy',ind);
atomdata.omegax = readposition(var,'omegax',ind);
atomdata.omegay = readposition(var,'omegay',ind);
atomdata.omegaz = readposition(var,'omegaz',ind);
atomdata.fx = readposition(var,'fx',ind);
atomdata.fy = readposition(var,'fy',ind);
atomdata.tqx = readposition(var,'tqx',ind);
atomdata.tqy = readposition(var,'tqy',ind);
atomdata.tqz = readposition(var,'tqz',ind);
atomdata.c2b1 = readposition(var,'c_2b[1]',ind);
atomdata.c2b2 = readposition(var,'c_2b[2]',ind);
atomdata.c2b4 = readposition(var,'c_2b[4]',ind);
atomdata.c2a1 = readposition(var,'c_2a[1]',ind);
atomdata.c2a2 = readposition(var,'c_2a[2]',ind);
atomdata.c2a3 = readposition(var,'c_2a[3]',ind);
atomdata.c2a4 = readposition(var,'c_2a[4]',ind);

function [k,ok] = findposition(structvar,varname)
    k = 1; ok = 0;
    while k <= length(structvar.header{1})
       if strcmp(structvar.header{1}{k},varname)
          ok = 1;
          break;
       else
          k = k+1;
       end
    end
end
function data = readposition(structvar,varname,in)
    k = 1; ok = 0;
    while k <= length(structvar.header{1})
       if strcmp(structvar.header{1}{k},varname)
          ok = 1;
          break;
       else
          k = k+1;
       end
    end
    if ok==1
        data = structvar.data(in,k);
    else
        data = zeros(size(atomdata.r));
    end
end
end
