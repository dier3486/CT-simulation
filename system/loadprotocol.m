function SYS = loadprotocol(SYS)
% load the protocol (of current series) to system

% protocol
protocol = SYS.protocol;
% I know the protocol is configure.protocol.series{ii}

% initial
1;

% tube
tube_corr = SYS.source.tube_corr;
% focal position
focalposition = reshape(tube_corr.focalposition, [], 3);
if isfield(tube_corr, 'focaldistort')
    fdis = tube_corr.focaldistort; 
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
SYS.source.focalposition = focalposition(protocol.focalspot, :);
SYS.source.focalnumber = size(SYS.source.focalposition, 1);
% KV mA (multi)
N_KV = length(protocol.KV);
N_mA = length(protocol.mA);
N_mAair = length(protocol.mA_air);
SYS.source.Wnumber = max(N_KV, N_mA);
if N_KV>1
    SYS.source.KV = num2cell(protocol.KV);
else
    SYS.source.KV = cell(SYS.source.Wnumber, 1);
    SYS.source.KV(:) = {protocol.KV};
end
if N_mA>1
    SYS.source.mA = num2cell(protocol.mA);
else
    SYS.source.mA = cell(SYS.source.Wnumber, 1);
    SYS.source.mA(:) = {protocol.mA};
end
if N_mAair>1
    SYS.source.mA_air = num2cell(protocol.mA_air);
else
    SYS.source.mA_air = cell(SYS.source.Wnumber, 1);
    SYS.source.mA_air(:) = {protocol.mA_air};
end
% focal size
focalsize = reshape(tube_corr.focalsize, [], 2);
SYS.source.focalsize = focalsize(protocol.focalsize, :);
% spectrum
SYS.source.spectrum = cell(SYS.source.Wnumber, 1);
tube_corr.main = reshape(tube_corr.main, [], tube_corr.KVnumber);
if strcmpi(SYS.simulation.spectrum, 'Single')
    samplekeV = SYS.world.referencekeV;
else
    samplekeV = SYS.world.samplekeV;
end
for ii = 1:SYS.source.Wnumber
    KV_ii = SYS.source.KV{ii};
    spectdata = reshape(tube_corr.main(:, (tube_corr.KVtag == KV_ii)), [], 2);
    SYS.source.spectrum{ii} = interp1(spectdata(:,1), spectdata(:,2), samplekeV);
    SYS.source.spectrum{ii}(isnan(SYS.source.spectrum{ii})) = 0;
end
% other
if isfield(tube_corr, 'offfocalintensity')
    SYS.source.offfocalintensity = tube_corr.offfocalintensity;
    SYS.source.offfocalwidth = tube_corr.offfocalwidth;
end
    

% collimation
% bowtie
% I know the bowtie index is
switch lower(protocol.bowtie)
    case {'empty', 0}
        % do nothing
        bowtie_index = [];
    case {'body', 'large', 1}
        bowtie_index = [1 2];
    case {'head', 'small', 2}
        bowtie_index = [1 3];
    otherwise
        error(['Unknown bowtie: ' protocol.bowtie]);
end
% loop multi-bowtie
for ii = 1:length(SYS.collimation.bowtie(:))
    % bowtie curve
    SYS.collimation.bowtie{ii} = ...
        getbowtiecurve(SYS.collimation.bowtie{ii}, SYS.source, bowtie_index);
    % bowtie material
    SYS.collimation.bowtie{ii}.material = SYS.collimation.bowtie{ii}.bowtie_corr.material;
end

% detector
% collimator -> detector 
if isfield(SYS, 'console') && isfield(SYS.console.protocaltrans, 'collimatorexplain')
    collimatorexplain = SYS.console.protocaltrans.collimatorexplain;
    % I know the collimatorexplain shall contian a cell field 'collimator'
else
    collimatorexplain = [];
end
SYS.detector = collimatorexposure(protocol.collimator, SYS.detector, SYS.detector.detector_corr, collimatorexplain);
% extra detector info
SYS.detector = detectorextracollim(SYS.detector, SYS.detector.detector_corr, samplekeV);

% DCB
SYS.datacollector.integrationtime = protocol.integrationtime;

% output
% reset output path
if isfield(protocol, 'outputpath')
    SYS.output.path = protocol.outputpath;
end
% output file names and version
SYS.output = outputfilenames(SYS.output, protocol, SYS.source);
% output style
if isfield(protocol, 'rawdatastyle') && ~isempty(protocol.rawdatastyle)
    SYS.output.rawdatastyle = protocol.rawdatastyle;
else
    % default style 
    SYS.output.rawdatastyle = '24bit';
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
% negative
% I know
sign_bowtiebox = sign(bowtie_corr.box(3));
% which should be sign(det(diag(bowtie_corr.box))), anyhow;
bowtie.bowtiecurve = bowtie.bowtiecurve.*sign_bowtiebox;

end


function output = outputfilenames(output, protocol, source)

files = struct();
corrversion = struct();
% namekey
if isfield(protocol, 'namekey') && ~isempty(protocol.namekey)
    namekey = ['_' protocol.namekey];
elseif isfield(output, 'namekey') && ~isempty(output.namekey)
    namekey = ['_' output.namekey];
else
    namekey = '';
end

% ini
Nw = source.Wnumber;
files.rawdata = cell(1, Nw);

% corr tables to output
if isfield(output, 'corrtable') && ~isempty(output.corrtable)
    corrtables = regexp(output.corrtable, '(, +)|(,)', 'split');
    corrtables = regexp(corrtables, '_', 'split');
    % ini corr
    for icorr = 1:length(corrtables)
        crr_i = corrtables{icorr};
        files.(crr_i{1}) = cell(1, Nw);
        if length(crr_i) > 1
            corrversion.(crr_i{1}) = crr_i{2};
        else
            corrversion.(crr_i{1}) = 'v1.0';
        end
    end
else
    corrtables = {};
end
% files to output
for iw = 1:Nw
    % name tags
    nametags = nametagrule(output.namerule, protocol, source.KV{iw}, source.mA{iw});
    % rawdata
    files.rawdata{iw} = ['rawdata' namekey nametags '_' output.rawdataversion];
    % corr table
    for icorr = 1:length(corrtables)
        table = corrtables{icorr}{1};
        files.(table){iw} = [table namekey nametags '_' corrversion.(table)];
    end
end
% NOTE: those names without EXT, .e.g. .raw or .corr.

% recon xml
files.reconxml = ['recon' namekey nametagrule(output.namerule, protocol)];

% to return
output.files = files;
output.corrversion = corrversion;
end
