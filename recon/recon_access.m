function [dataflow, prmflow, status] = recon_access(status, echo_onoff, dataflow, prmflow)
% recon & cali governing function
% [dataflow, prmflow, status] = recon_access(status)

if nargin<2
    echo_onoff = false;
end
if nargin<3
    % initial
    dataflow = struct();
    prmflow = struct();
end

% initial steps
if echo_onoff, fprintf('Recon Series %d\n', status.series_index); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');
if ~status.jobdone
    return
end

if echo_onoff, fprintf('  load calibration tables...'); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadcorrs');
if ~status.jobdone
    if echo_onoff, fprintf(' failed\n'); end
    return
else
    if echo_onoff, fprintf(' done\n'); end
end

% load rawdata
if echo_onoff, fprintf('  read rawdata...'); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadrawdata');
if ~status.jobdone
    if echo_onoff, fprintf(' failed\n'); end
    return
else
    if echo_onoff, fprintf(' done\n'); end
end
% for large data we should employ view buffer in loading data, TBC

% run pipe nodes
pipefields = fieldnames(prmflow.pipe);
for i_node = 1:length(pipefields)
    node = pipefields{i_node};
    if echo_onoff, fprintf('  [recon node] %s...', node); end
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, node);
    if ~status.jobdone
        if echo_onoff, fprintf(' failed\n'); end
        return
    else
        if echo_onoff, fprintf(' done\n'); end
    end
end
if echo_onoff, fprintf('Done\n'); end

end
