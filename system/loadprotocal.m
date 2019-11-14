function SYS = loadprotocal(SYS)
% load the protocal (of current series) to system

% protocal
protocal = SYS.protocal;
% I know the protocal is configure.protocal.series{ii}

% initial

% detector
% prepare detector_corr
det_corr = SYS.detector.detector_corr;
det_corr.position = reshape(det_corr.position, det_corr.Npixel, det_corr.Nslice, 3);
% values to copy
SYS.detector.Npixel = det_corr.Npixel;
% collimator -> detector 
SYS.detector = collimatorondetector(protocal.collimator, SYS.detector, det_corr);

% tube
% focal position
focalposition = reshape(SYS.source.tube_corr.focalposition, [], 3);
if isfield(SYS.source.tube_corr, 'focaldistort')
    fdis = SYS.source.tube_corr.focaldistort; 
    if length(fdis(:)) <= 1
        focalposition(1) = focalposition(1) + fdis;
    else
        fdis = reshape(fdis, [], 3);
        if size(focalposition, 1) == 1
            focalposition = repmat(focalposition, size(fdis,1), 1);
        end
        focalposition = focalposition + fdis;
    end
end
SYS.source.focalposition = focalposition(protocal.focalspot, :);
SYS.source.focalnumber = size(SYS.source.focalposition, 1);
% focal size
focalsize = reshape(SYS.source.tube_corr.focalsize, [], 2);
SYS.source.focalsize = focalsize(protocal.focalsize, :);
% KV mA (multi)
N_KV = length(protocal.KV);
N_mA = length(protocal.mA);
SYS.source.Wnumber = max(N_KV, N_mA);
if N_KV>1
    SYS.source.KV = num2cell(protocal.KV);
else
    SYS.source.KV = cell(SYS.source.Wnumber, 1);
    SYS.source.KV(:) = {protocal.KV};
end
if N_mA>1
    SYS.source.mA = num2cell(protocal.mA);
else
    SYS.source.mA = cell(SYS.source.Wnumber, 1);
    SYS.source.mA(:) = {protocal.mA};
end
% spectrum
SYS.source.spectrum = cell(SYS.source.Wnumber, 1);
tube_corr = SYS.source.tube_corr;
tube_corr.main = reshape(tube_corr.main, [], tube_corr.KVnumber);
samplekeV = SYS.world.samplekeV;
for ii = 1:SYS.source.Wnumber
    KV_ii = SYS.source.KV{ii};
    spectdata = reshape(tube_corr.main(:, (tube_corr.KVtag == KV_ii)), [], 2);
    SYS.source.spectrum{ii} = interp1(spectdata(:,1), spectdata(:,2), samplekeV);
    SYS.source.spectrum{ii}(isnan(SYS.source.spectrum{ii})) = 0;
end

% collimator
% bowtie
switch lower(protocal.bowtie)
    case 'empty'
        % do nothing
        bowtie_index = [];
    case 'body'
        bowtie_index = [1 2];
    case 'head'
        bowtie_index = [3 4];
    otherwise
        error(['Unknown bowtie: ' protocal.bowtie]);
end
SYS.collimation.bowtie = ...
    getbowtiecurve(SYS.collimation.bowtie, SYS.source, bowtie_index);

% output
% output file names
SYS.output.files = outputfilenames(SYS.output, protocal);

end


function detector = collimatorondetector(collimator, detector, det_corr)
% explain th ecollimator (hard code)

switch collimator
    case {'16x0.55', '16x0.5'}
        sliceindex = 5:20;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 16;
        detector.startslice = 5;
        detector.endslice = 20;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:16;
    case {'16x1.1', '16x1.0'}
        detector.position = reshape(det_corr.position, [], 3);
        detector.Nslice = 24;
        detector.startslice = 1;
        detector.endslice = 24;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO*2;
        detector.slicemerge = [1 1 2 2 3 3 4 4 5 6 7 8 9 10 11 12 13 13 14 14 15 15 16 16];
    case {'8x1.1', '8x1.0'}
        sliceindex = 5:20;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 16;
        detector.startslice = 5;
        detector.endslice = 20;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO*2;
        detector.slicemerge = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
    case {'8x0.55', '8x0.5'}
        sliceindex = 8:15;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 8;
        detector.startslice = 8;
        detector.endslice = 15;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:8;
    case {'4x0.55', '4x0.5'}
        sliceindex = 10:13;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 4;
        detector.startslice = 10;
        detector.endslice = 13;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:4;
    otherwise
        error(['Illeagal collimator: ', protocal.collimator]);  
end
end


function bowtie = getbowtiecurve(bowtie, source, bowtie_index)
% get bowtie curve

if nargin < 3
    bowtie_index = [];
end

if isempty(bowtie_index)
    % empty bowtie
    bowtie.anglesample = [];
    bowtie.bowtiecurve = [];
    return
end

% original focal position
focal_orig = source.tube_corr.focalposition(1,:);
% bowtie_corr
bowtie_corr = bowtie.bowtie_corr;
% use bowtie_index to select bowtie
bowtie_corr.main = reshape(bowtie_corr.main, bowtie_corr.Nsample, []);
bowtiecv_orig = bowtie_corr.main(:, bowtie_index);
% I know the axis-x is from 0 to bowtie_corr.box(1)
bowtiecv_orig(:, 1) = bowtiecv_orig(:, 1) - bowtie_corr.box(1)/2;
% initial the returns
bowtie.bowtiecurve = zeros(bowtie_corr.Nsample, source.focalnumber);
bowtie.anglesample = zeros(bowtie_corr.Nsample, source.focalnumber);
% to loop the focal spots
for i_focal = 1:source.focalnumber
    focaltobottom = bowtie_corr.focaltobottom + ...
        focal_orig(1)- source.focalposition(i_focal, 1);
    bowtie.anglesample(:, i_focal) = atan2(focaltobottom - bowtiecv_orig(:, 2), ...
        bowtiecv_orig(:, 1) - source.focalposition(i_focal, 1)) - pi/2;
    bowtie.bowtiecurve(:, i_focal) = bowtiecv_orig(:, 2).*sec(bowtie.anglesample(:, i_focal));
end

end


function files = outputfilenames(output, protocal)

files = struct();
% namekey
namekey = output.namekey;
if ~isempty(namekey)
    namekey = ['_' namekey];
end

switch lower(output.namerule)
    case {'default'}
        % rawdata
        rawtags = ['_series' num2str(protocal.series_index) '_' ...
            protocal.scan '_' protocal.bowtie '_' protocal.collimator ...
            '_' num2str(protocal.KV) 'KV' num2str(protocal.mA) 'mA' '_' ...
            num2str(protocal.rotationspeed) 'secprot'];
        files.rawdata = ['rawdata' namekey rawtags];
        % air 
        if strfind(output.corrtable, 'air')
            airtags = ['_' protocal.bowtie '_' protocal.collimator '_' ...
                num2str(protocal.KV) 'KV' num2str(protocal.mA_air) 'mA'];
            files.aircorr = ['air' airtags];
        end
        % offset
        % TBC
    otherwise
        % most simple filenames
        rawtags = ['_series' num2str(protocal.series_index)];
        files.rawdata = ['rawdata' namekey rawtags];
        % air
        if strfind(output.corrtable, 'air')
            files.aircorr = 'air';
        end
        % offset
        % TBC
end
% NOTE: those names without version tag, e.g. _v1.0, and EXT, .e.g. .raw.
end