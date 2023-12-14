function [datastruct, cfg] = loadbindata(datafile, cfgfile, warnonoff)
% load bin data file

if nargin<3
    warnonoff = true;
end

if isstruct(cfgfile)
    cfg = cfgfile;
elseif isfile(cfgfile)
    cfg = readcfgfile(cfgfile);
elseif exist(cfgfile, 'file')
    if warnonoff
        warning('To load configure file in MATLAB search path! It might be unsafe.');
    end
    cfg = readcfgfile(cfgfile);
else
    error('Can not open the cfgfile!');
end

fid = fopen(datafile, 'r');

datastruct = sparsepack(fid, cfg);
fclose(fid);

end
