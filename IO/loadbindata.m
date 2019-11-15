function [datastruct, cfg] = loadbindata(datafile, cfgfile)
% load bin data file

cfg = readcfgfile(cfgfile);

fid = fopen(datafile, 'r');
data = fread(fid, inf, 'uint8=>uint8');
fclose(fid);
datastruct = sparsepack(data, cfg);

end
