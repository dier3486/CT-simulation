function [dataflow, prmflow, status] = recon_access(status, echo_onoff)
% recon & cali governing function
% [dataflow, prmflow, status] = recon_access(status)

if nargin<2
    echo_onoff = false;
end

% initial
dataflow = struct();
prmflow = struct();
% initial steps
if echo_onoff, fprintf('Recon Series %d\n', status.series_index); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');
if echo_onoff, fprintf('  load calibration tables...'); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadcorrs');
if echo_onoff, fprintf(' done\n'); end

% load rawdata
if echo_onoff, fprintf('  read rawdata...'); end
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadrawdata');
if echo_onoff, fprintf(' done\n'); end
% for large data we should employ view buffer in loading data, TBC

% run pipe nodes
pipefields = fieldnames(prmflow.pipe);
for i_node = 1:length(pipefields)
    node = pipefields{i_node};
    if echo_onoff, fprintf('  [recon node] %s...', node); end
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, node);
    if echo_onoff, fprintf(' done\n'); end
end
if echo_onoff, fprintf('Done\n'); end

end
