function prm = parameterforcorr(SYS, corrversion)
% return the paramters for simu corr tables from SYS and SYS.protocol in simulation
% prm = parameterforcorr(SYS, corrversion);

% NOTE: plz self-discipline this code as same as caliprmforcorr.m

if nargin<2
    corrversion = 'v1.0';
end

% ID
versionID = regexp(corrversion, '\d+','match');
ID = [0 0 0 0];
for ii = 1:min(length(versionID), 4)
    ID(end-ii+1) = str2double(versionID{end-ii+1});
end
prm.ID = ID;

% Series_Number
prm.seriesnumber = SYS.protocol.series_index;

% Npixel
prm.Npixel = SYS.detector.Npixel;

% Start Slice
prm.startslice = SYS.detector.startslice;

% End Slice
prm.endslice = SYS.detector.endslice;

% Slice merge
prm.slicemerge = SYS.detector.mergescale;

% Slice Number
prm.slicenumber = max(SYS.detector.slicemerge);
% NOTE: which is not SYS.detector.Nslice

% view number
prm.viewnumber = SYS.protocol.viewnumber;

% focalspot
focalspot = SYS.protocol.focalspot;
prm.focalspot = sum(2.^(focalspot-1));
% I know the rule is that

% focalsize
prm.focalsize = SYS.protocol.focalsize;

% bowtie
bowtie = SYS.protocol.bowtie;
if isnumeric(bowtie)
    prm.bowtie = bowtie;
else
    switch lower(bowtie)
        case 'empty'
            prm.bowtie = 0;
        case {'body', 'large'}
            prm.bowtie = 1;
        case {'head', 'small'}
            prm.bowtie = 2;
        otherwise
            % unknown bowtie
            prm.bowtie = -1;
    end
end

% KV
prm.KV = SYS.source.KV;

% mA
prm.mA = SYS.source.mA;

% mA_air
prm.mA_air = SYS.source.mA_air;

% focalnumber
prm.focalnumber = SYS.source.focalnumber;

% rotationspeed
prm.rotationspeed = SYS.protocol.rotationspeed;

end

