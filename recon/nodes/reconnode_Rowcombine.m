function [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status)
% recon node, row combine
% [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status);

% parameters to use in prmflow
% Nview = prmflow.recon.Nview;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;

if isfield(prmflow.system.detector, 'rowcombine') && ~isempty(prmflow.system.detector, 'rowcombine')
    rowcombine = prmflow.system.detector.rowcombine;
    % combine the rawdata
    [dataflow.rawdata, Nrowcombined] = detectorslicemerge(dataflow.rawdata, Npixel, Nslice, rowcombine, 'mean');
    % combine the detector position
    [prmflow.system.detector.position, ~] = ...
        detectorslicemerge(prmflow.system.detector.position, Npixel, Nslice, rowcombine, 'mean');
    prmflow.recon.Npixel = Nrowcombined;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end