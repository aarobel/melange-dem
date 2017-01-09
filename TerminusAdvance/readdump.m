function results = readdump(varargin)
% readdump -- reads data from a LAMMPS dump file. 
%
% Author: Agnieszka Herman, IOUG (agnieszka.herman@ug.edu.pl)
%
dump = fopen(varargin{1},'r');
if dump<=0
    error('Dumpfile not found!');
end

n = 1;
while feof(dump) == 0 
    id = fgetl(dump);
    if strncmpi(id,'ITEM: TIMESTEP',14)
        t(n) = str2num(fgetl(dump));                            %#ok<*ST2NM>
    elseif strncmpi(id,'ITEM: NUMBER OF ATOMS',21)
        N(n) = str2num(fgetl(dump));
    elseif strncmpi(id,'ITEM: BOX BOUNDS',16)
        xb(n,:) = str2num(fgetl(dump));
        yb(n,:) = str2num(fgetl(dump));
        zb(n,:) = str2num(fgetl(dump));
    elseif strncmpi(id,'ITEM: ATOMS',11)
        tmp = strtrim(id(13:end));
        ib = [0 strfind(tmp,' ') length(tmp)+1];
        for k = 1:length(ib)-1
            header{n}{k} = tmp(ib(k)+1:ib(k+1)-1);
        end
        ncols = length(strfind(id,' ')) - 2;
        data(n,:,:) = fscanf(dump,'%f',[ncols N(n)])';        
        n = n+1;
    elseif strncmpi(id,'ITEM: NUMBER OF ENTRIES',23)
        Np(n) = str2num(fgetl(dump));
    elseif strncmpi(id,'ITEM: ENTRIES',13)
        tmp = strtrim(id(15:end));
        ib = [0 strfind(tmp,' ') length(tmp)+1];
        for k = 1:length(ib)-1
            headerentries{n}{k} = tmp(ib(k)+1:ib(k+1)-1);
        end
        ncols = length(strfind(id,' ')) - 2;
        entries{n} = fscanf(dump,'%f',[ncols Np(n)])';
        n = n+1;
    end
end
if exist('data','var')
    data = squeeze(data);
    results = struct('t',t','N',N','xb',xb,'yb',yb,'zb',zb,'data',data);
    results.header = header;
elseif exist('entries','var')
    entries = squeeze(entries);
    headerentries = squeeze(headerentries);
    results = struct('t',t','Np',Np');
    results.headerentries = headerentries;
    results.entries = entries;
end
fclose(dump);
