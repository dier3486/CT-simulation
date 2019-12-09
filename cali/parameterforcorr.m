function prm = parameterforcorr(SYS)
% return the paramters for corr tables form SYS and SYS.protocol
% prm = parameterforcorr(SYS)

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


