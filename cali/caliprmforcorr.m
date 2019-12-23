function prm = caliprmforcorr(prmflow, corrversion)
% return the paramters for corr tables from protocol in calibration
% prm = caliprmforcorr(protocol, corrversion)

if nargin<2
    corrversion = 'v1.0';
end
% to use
system = prmflow.system;
protocol = prmflow.protocol;

% ID
versionID = regexp(corrversion, '\d+','match');
ID = [0 0 0 0];
for ii = 1:min(length(versionID), 4)
    ID(end-ii+1) = str2double(versionID{end-ii+1});
end
prm.ID = ID;

% Npixel
prm.Npixel = system.detector.Npixel;

% Start Slice
prm.startslice = system.detector.startslice;

% End Slice
prm.endslice = system.detector.endslice;

% Slice merge
prm.slicemerge = system.detector.mergescale;

% Slice Number
prm.slicenumber = system.detector.Npixel;

% view number
prm.viewnumber = protocol.viewnumber;

% focalspot
prm.focalspot = protocol.focalspot;

% bowtie
bowtie = protocol.bowtie;
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
prm.KV = protocol.KV;

% mA
prm.mA = protocol.mA;

% rotationspeed
prm.rotationspeed = protocol.rotationspeed;

% focalnumber
prm.focalnumber = system.Nfocal;

end

