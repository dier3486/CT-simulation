function [prmflow, status] = reconinitial(prmflow, status)
% recon initial
% [prmflow, status] = reconinitial(prmflow, status)

% copy status.reconcfg to prmflow
if ~iscell(status.reconcfg)
    reconcfg = status.reconcfg;
else
    reconcfg = status.reconcfg{status.series_index};
end

if isempty(reconcfg)
    % empty configure?
    status.jobdone = false;
    status.errorcode = 1;
    status.errormsg = '[reconinitial] empty recon configure';
    return;
end

% copy reconcfg to prmflow
prmflow = structmerge(reconcfg, prmflow, 0, 0);
% but maintain the extra fields in prmflow

% reload sub-config file
prmflow = subconfigure(prmflow);

% external supports
if isfield(prmflow, 'external')
    prmflow = CRIS2prmflow(prmflow, prmflow.external.rawxml);
end

% clean
prmflow = iniprmclean(prmflow);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function prmflow = iniprmclean(prmflow)
% to fill up the paramters which could be used in recon but not configured
% hard code

% collimatorexplain
if ~isfield(prmflow.system, 'collimatorexplain')
    prmflow.system.collimatorexplain = [];
elseif ischar(prmflow.system.collimatorexplain)
    prmflow.system.collimatorexplain = readcfgfile(prmflow.system.collimatorexplain);
end
% IOstandard
if ~isfield(prmflow, 'IOstandard')
    prmflow.IOstandard = [];
end
% ini corrtable
if ~isfield(prmflow, 'corrtable')
    prmflow.corrtable = struct();
end
% ini recon (always)
prmflow.recon = struct();

% explain focal spot
spots = fliplr(dec2bin(prmflow.protocol.focalspot)=='1');
prmflow.system.Nfocal = sum(spots);
prmflow.system.focalspot = find(spots);

end