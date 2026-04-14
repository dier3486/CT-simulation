function [viewnum, viewsize] = checksize(datafile, cfgpath)
% check the view number of a rawdatafile

if nargin<2
    cfgpath = '';
    warnonoff = false;
else
    warnonoff = true;
end

% to char
datafile= char(datafile);

% file format configure
if isfile(cfgpath)
    cfgfile = cfgpath;
else
    cfgfile = cfgmatchrule(datafile, cfgpath); 
end
if isfile(cfgfile)
    cfg = readcfgfile(cfgfile);
elseif exist(cfgfile, 'file')
    if warnonoff
        warning('To load configure file in MATLAB search path! It might be unsafe.');
    end
    cfg = readcfgfile(cfgfile);
else
    error('Can not open the cfgfile!');
end

% open the file
fid = fopen(datafile, 'r');

% check datasize
datasize = fileposorgtoend(fid);

% decode view size
cfg.size = decodenumber(cfg.size, []);
cfg.number = 1;

if ~isavail(cfg.size)
    % size is unknown
    % try to read first view    
    [~, cfg] = sparsepack(fid, cfg);
end

% fclose
fclose(fid);

% to return
viewnum = floor(datasize/cfg.size);
viewsize = cfg.size;

end