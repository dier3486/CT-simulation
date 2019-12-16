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
    detectorslicemerge(prmflow.system.detector.position, prmflow.system.detector, 'mean');
prmflow.system.detector.Nmergedslice = Nmergedslice;
% copy other parameters from det_corr
prmflow.system.detector = structmerge(prmflow.system.detector, det_corr);

% to prm.recon
prmflow.recon.Nslice = Nmergedslice;
prmflow.recon.Npixel = double(prmflow.system.detector.Npixel);

% other tables
pipenodes = fieldnames(prmflow.pipe);
for ii = 1:length(pipenodes)
    if isfield(prmflow.pipe.(pipenodes{ii}), 'corr')
        prmflow.corrtable.(pipenodes{ii}) = loaddata(prmflow.pipe.(pipenodes{ii}).corr, prmflow.IOstandard);
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


