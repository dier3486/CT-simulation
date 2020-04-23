function [prmflow, status] = loadcalitables(prmflow, status)
% load calibration tables
% [prmflow, status] = loadcalitables(reconcfg, prmflow, status);

% detector position
% load corr file
detcorrfile = prmflow.system.detector_corr;
det_corr = loaddata(detcorrfile, prmflow.IOstandard);
% explain the collimator
prmflow.system.detector = collimatorexposure(prmflow.protocol.collimator, ...
    [], det_corr, prmflow.system.collimatorexplain);
% mergeslice
prmflow.system.detector.position = reshape(prmflow.system.detector.position, [], 3);
[prmflow.system.detector.position, Nmergedslice] = ...
    detectorslicemerge(prmflow.system.detector.position, prmflow.system.detector.Npixel, prmflow.system.detector.Nslice, ...
    prmflow.system.detector.slicemerge, 'mean');
prmflow.system.detector.Nmergedslice = Nmergedslice;
% copy other parameters from det_corr
prmflow.system.detector = structmerge(prmflow.system.detector, det_corr);

% to prm.recon
prmflow.recon.Nslice = Nmergedslice;
prmflow.recon.Npixel = double(prmflow.system.detector.Npixel);

% other tables
% corrpath is the path to looking for corr files
if isfield(prmflow.system, 'corrpath') && ~isempty(prmflow.system.corrpath)
    corrpath = prmflow.system.corrpath;
else
    [corrpath, ~, ~] = fileparts(prmflow.rawdata);
end
% corrext is the ext of the corr files
if isfield(prmflow.system, 'corrext')
    corrext = prmflow.system.corrext;
    % NOTE: set corrext = '.(corr|ct)' to looking for both .corr and .ct files
else
    % default corrext is .corr
    corrext = '.corr';
end
pipenodes = fieldnames(prmflow.pipe);
for ii = 1:length(pipenodes)
    if isfield(prmflow.pipe.(pipenodes{ii}), 'corr')
        % if corr is empty, auto looking for a mtached corr file
        if isempty(prmflow.pipe.(pipenodes{ii}).corr)
            % name of the corr
            corrname = regexp(pipenodes{ii}, '_', 'split');
            corrname = lower(corrname{1});
            % looking for corr file
            corrfile = corrcouplerule(prmflow.protocol, corrpath, prmflow.system.filematchrule, corrname, corrext);
            prmflow.pipe.(pipenodes{ii}).corr = corrfile;
        end
        % load the corrfile
        prmflow.corrtable.(pipenodes{ii}) = loaddata(prmflow.pipe.(pipenodes{ii}).corr, prmflow.IOstandard);
    end
    % else, no corr to load for this node
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


