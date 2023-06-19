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
if isfield(tube_corr, 'tubenumber')
    % multi-tube CT
    SYS.source.tubenumber = tube_corr.tubenumber;
    SYS.source.origfocalpos = focalposition;
    focaloffset = reshape(tube_corr.focaloffset, [], 3);
    focaloffset = focaloffset(protocol.focalspot, :);
    SYS.source.focalposition = reshape(reshape(focalposition, 1, [], 3) +reshape(focaloffset,[], 1, 3), [], 3);
    SYS.source.focalnumber = size(focaloffset, 1);
else
    % single-tube CT
    SYS.source.tubenumber = 1;
    SYS.source.origfocalpos = focalposition(1, :);
    SYS.source.focalposition = focalposition(protocol.focalspot, :);
    SYS.source.focalnumber = size(SYS.source.focalposition, 1);
end

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
    SYS.source.spectrum{ii} = spectrumresample(spectdata(:,1), spectdata(:,2), samplekeV);
    SYS.source.spectrum{ii}(isnan(SYS.source.spectrum{ii})) = 0;
end
% offfocal
if isfield(tube_corr, 'offfocal')
    if isfield(SYS.source, 'offfocal')
        SYS.source.offfocal = structmerge(SYS.source.offfocal, tube_corr.offfocal);
    else
        SYS.source.offfocal = tube_corr.offfocal;
    end
end

% collimation
% bowtie
% I know the bowtie index is
switch lower(protocol.bowtie)
    case {'empty', 'air', 0}
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
if isfield(SYS.collimation, 'bowtie')
    for ii = 1:length(SYS.collimation.bowtie(:))
        % bowtie curve
        SYS.collimation.bowtie{ii} = ...
            getbowtiecurve(SYS.collimation.bowtie{ii}, SYS.source, bowtie_index);
        % bowtie material
        SYS.collimation.bowtie{ii}.material = SYS.collimation.bowtie{ii}.bowtie_corr.material;
    end
end
% filter
if isfield(SYS.collimation, 'filter')
    for ii = 1:length(SYS.collimation.filter(:))
        % original focal position
        origfocalpos = SYS.source.origfocalpos(1,:);
        focalpos = SYS.source.focalposition(1:SYS.source.focalnumber, :);
        SYS.collimation.filter{ii}.origangle = atan2(-focalpos(:,2), -focalpos(:,1)) - ...
            atan2(-origfocalpos(2), -origfocalpos(1));
    end
end

% detector
% collimator -> detector 
if isfield(SYS, 'console') && isfield(SYS.console.protocoltrans, 'collimatorexplain') && ...
        ~isempty(SYS.console.protocoltrans.collimatorexplain)
    collimatorexplain = SYS.console.protocoltrans.collimatorexplain.collimator;
    % I know the collimatorexplain shall contian a cell field 'collimator'
else
    collimatorexplain = [];
end
SYS.detector = collimatorexposure(protocol.collimator, SYS.detector, SYS.detector.detector_corr, collimatorexplain);
% extra detector info
SYS.detector = detectorextracollim(SYS.detector, SYS.detector.detector_corr, samplekeV);
% copy other parameters from det_corr
SYS.detector = structmerge(SYS.detector, SYS.detector.detector_corr);
% pixel range
if isfield(SYS.detector, 'pixelrange')
    SYS.detector.pixelrange = reshape(SYS.detector.pixelrange, 2, []);
    SYS.detector.Nprange = max(mod(SYS.detector.pixelrange(2, :)-SYS.detector.pixelrange(1, :), SYS.detector.Npixel) + 1);
end

% DCB
SYS.datacollector.integrationtime = protocol.integrationtime;
% log2?
if isfield(protocol, 'rawdatastyle')
    switch protocol.rawdatastyle
        case '16bit'
            SYS.datacollector.islog2 = true;
        case '24bit'
            SYS.datacollector.islog2 = false;
        otherwise
            % default value was set in systemconfigure.m
            1;
    end
end

% reset output path
if isfield(protocol, 'outputpath')
    SYS.output.path = protocol.outputpath;
end

% console
% name rule
if isfield(SYS.console.protocoltrans, 'filetagsrule')
    filetagsrule = SYS.console.protocoltrans.filetagsrule;
else
    filetagsrule = [];
end
% output file names and version
SYS.output = outputfilenames(SYS.output, protocol, SYS.source, filetagsrule);
% output style
if isfield(protocol, 'rawdatastyle') && ~isempty(protocol.rawdatastyle)
    SYS.output.rawdatastyle = protocol.rawdatastyle;
else
    % default style
    if SYS.datacollector.islog2
        SYS.output.rawdatastyle = '16bit';
    else
        SYS.output.rawdatastyle = '24bit';
    end
end
% filematchrule
% nothing to do



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
if isfield(source, 'origfocalpos')
    origfocalpos = source.origfocalpos(1,:);
else
    origfocalpos = source.tube_corr.focalposition(1,:);
end
% origfocalang = atan2(-origfocalpos(2), -origfocalpos(1));
% bowtie_corr
bowtie_corr = bowtie.bowtie_corr;
% bowtie distort
if isfield(bowtie_corr, 'bowtiedistort') && ~isempty(bowtie_corr.bowtiedistort)
    bowtiedistort = bowtie_corr.bowtiedistort;
else
    bowtiedistort = [0 0 0];
end
% bowtiedistort = [0 0 0]; 
% use bowtie_index to select bowtie
bowtie_corr.main = reshape(bowtie_corr.main, bowtie_corr.Nsample, []);
bowtiecv_orig = bowtie_corr.main(:, bowtie_index);
% I know the axis-x is from 0 to bowtie_corr.box(1)
bowtiecv_orig(:, 1) = bowtiecv_orig(:, 1) - bowtie_corr.box(1)/2;
% % debug
% bowtiecv_orig(:, 2) = -bowtiecv_orig(:, 2);
% initial the returns
bowtie.bowtiecurve = zeros(bowtie_corr.Nsample, source.focalnumber);
bowtie.anglesample = zeros(bowtie_corr.Nsample, source.focalnumber);
% loop the focal spots
for i_focal = 1:source.focalnumber
    % to calculate the bowtie sample angles and bowtiecurve
    Dorig = sqrt(sum(origfocalpos(1:2).^2));
    Vnorigy = origfocalpos(1:2)'./Dorig;
    Vnorigx = [Vnorigy(2); -Vnorigy(1)];
    rfocal = source.focalposition(i_focal, 1:2)-origfocalpos(1:2);
    focaltobottom = bowtie_corr.focaltobottom + rfocal*Vnorigy;
    % I know the bowtie bottom is vertical to the line from ISO center to
    % the origfocalpos.
    y_bowtie = focaltobottom - bowtiecv_orig(:, 2) + bowtiedistort(2);
    x_bowtie = bowtiecv_orig(:, 1) + rfocal*Vnorigx + bowtiedistort(1);
    anlge_bowtie = atan2(y_bowtie, x_bowtie) - pi/2;
    % It could be a small angle between the focalposition to origfocalpos,
    % which reads,
    angleshift = atan2(-origfocalpos(2), -origfocalpos(1)) - ...
        atan2(-source.focalposition(i_focal,2), -source.focalposition(i_focal,1));
    % mod it to [-pi, pi).
    angleshift = mod(angleshift+pi, pi*2) - pi;
    % the anglesample of an X-ray path is the angle from the mid-ray to
    % that ray.
    bowtie.anglesample(:, i_focal) = anlge_bowtie + angleshift;
    bowtie.bowtiecurve(:, i_focal) = bowtiecv_orig(:, 2).*sec(anlge_bowtie);
end
% negative
% I know
sign_bowtiebox = sign(bowtie_corr.box(3));
% which should be sign(det(diag(bowtie_corr.box))), anyhow;
% % debug
% sign_bowtiebox = - sign_bowtiebox;
bowtie.bowtiecurve = bowtie.bowtiecurve.*sign_bowtiebox;

end


function output = outputfilenames(output, protocol, source, filetagsrule)

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
    if isfield(filetagsrule, 'rawdata')
        ruletags = filetagsrule.rawdata;
    else
        ruletags = [];
    end
    nametags = nametagrule(output.namerule, protocol, ruletags, source.KV{iw}, source.mA{iw});
    % rawdata
    files.rawdata{iw} = ['rawdata' namekey nametags '_' output.rawdataversion];
    % corr table
    for icorr = 1:length(corrtables)
        table = corrtables{icorr}{1};
        if isfield(filetagsrule, table)
            ruletags = filetagsrule.(table);
        else
            ruletags = [];
        end
        nametags = nametagrule(output.namerule, protocol, ruletags, source.KV{iw}, source.mA{iw});
        files.(table){iw} = [table namekey nametags '_' corrversion.(table)];
    end
end
% NOTE: those names without EXT, .e.g. .raw or .corr.

% recon xml
if isfield(filetagsrule, 'reconxml')
    ruletags = filetagsrule.reconxml;
else
    ruletags = [];
end
files.reconxml = ['recon' namekey nametagrule(output.namerule, protocol, ruletags)];

% to return
output.files = files;
output.corrversion = corrversion;
end
