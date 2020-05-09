function [prmflow, status] = reconinitial(prmflow, status)
% recon initial
% [prmflow, status] = reconinitial(prmflow, status)

% copy status.reconcfg to prmflow
if ~iscell(status.reconcfg)
    reconcfg = status.reconcfg;
else
    reconcfg = status.reconcfg{status.seriesindex};
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
else
    error('Illegal prmflow.system.collimatorexplain class: %s!', class(prmflow.system.collimatorexplain));
end
% filetagsrule
if ~isfield(prmflow.system, 'filetagsrule')
    prmflow.system.filetagsrule = struct();
elseif ischar(prmflow.system.filetagsrule)
    prmflow.system.filetagsrule = readcfgfile(prmflow.system.filetagsrule);
else
    error('Illegal prmflow.system.filetagsrule class: %s!', class(prmflow.system.filetagsrule));
end
% filematchrule
if ~isfield(prmflow.system, 'filematchrule')
    prmflow.system.filematchrule = struct();
elseif ischar(prmflow.system.filematchrule)
    prmflow.system.filematchrule = readcfgfile(prmflow.system.filematchrule);
else
    error('Illegal prmflow.system.filematchrule class: %s!', class(prmflow.system.filematchrule));
end

% IOstandard
if ~isfield(prmflow, 'IOstandard')
    prmflow.IOstandard = [];
end
% ini corrtable (always)
% if ~isfield(prmflow, 'corrtable')
%     prmflow.corrtable = struct();
% end
prmflow.corrtable = struct();
% ini recon (always)
prmflow.recon = struct();
% NOTE: sometimes we need to maintain data in prmflow for follow-up series, so we don't clean most of the informations in 
% prmflow when they already exist. But the prmflow.recon will always be cleaned.



end