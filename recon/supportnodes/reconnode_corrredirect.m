function [dataflow, prmflow, status] = reconnode_corrredirect(dataflow, prmflow, status)
% virus node, replace prmflow.corrtable.xxx by dataflow.xxxcorr
% [dataflow, prmflow, status] = reconnode_corrredirect(dataflow, prmflow, status);
% NOTE: 1. the node loadcalitables.m could overwrite the prmflow.corrtable.*,
%       2. the node reconinitial.m could overwrite any information if they were included in status.reconcfg,
%       3. other nodes also could overwrite the data you want to maintain,
% therefore if you want to send information by prmflow to next series, you have to avoid them be overwriten by those nodes in a
% normal recon pipe line.

% parameters set in pipe
credprm = prmflow.pipe.(status.nodename);
if isfield(credprm, 'corrext')
    corrext = credprm.corrext;
else
    corrext = 'corr';
end

if isfield(credprm, 'nodes')
    if iscell(credprm.nodes)
        for ii = 1:length(credprm.nodes)
            nodename_slip = regexp(credprm.nodes{ii}, '_', 'split');
            corrname = [lower(nodename_slip{1}) corrext];
            if isfield(dataflow, corrname)
                prmflow.corrtable.(credprm.nodes{ii}) = dataflow.(corrname);
            end
        end
    else
        nodename_slip = regexp(credprm.nodes, '_', 'split');
        corrname = [lower(nodename_slip{1}) corrext];
        if isfield(dataflow, corrname)
            prmflow.corrtable.(credprm.nodes) = dataflow.(corrname);
        end
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

