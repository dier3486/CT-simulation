function [dataflow, prmflow, status] = recon_access(dataflow, prmflow, status)
% recon & cali governing function
% [dataflow, prmflow, status] = recon_access(dataflow, prmflow, status)

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


if ~isfield(status, 'echo_onoff')
    status.echo_onoff = true;
end
if ~isfield(status, 'warning_onoff')
    status.warning_onoff = true;
end

% go
echo_onoff = status.echo_onoff;
if echo_onoff, fprintf('Recon Series %d\n', status.seriesindex); end

% initial steps (and initial GPU)
if echo_onoff, fprintf('  initial (GPU)...'); end
tic;
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');
timecost = toc;
[dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
if failedflag
    return;
end

% load calibration tables
if echo_onoff, fprintf('  load calibration tables...'); end
tic;
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadcorrs');
timecost = toc;
[dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
if failedflag
    return;
end

% pipeline nodes
pipelinenodes = fieldnames(status.pipeline);
if isempty(pipelinenodes)
    % ??
    return
end

% pipeline prepare
if prmflow.protocol.pipelineprepare
    if echo_onoff, fprintf('  pipeline prepare...'); end
    tic;
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'pipelineprepare');
    timecost = toc;
    [dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
    if failedflag
        return;
    end
end

% pipeline kick off
node = pipelinenodes{1};
while ~status.seriesdone
    % echo
    if echo_onoff, fprintf('  [recon node] %s...', node); end
    tic;
    % node prepare
    if ~prmflow.protocol.pipelineprepare
        1;
    end
    % node main
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, node);
    timecost = toc;
    1;
    [dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
    if failedflag
        return;
    end

    nextnode = rollpoling(status.pipeline);
    status.seriesdone = isempty(nextnode);
    if strcmp(node, nextnode)
        pause(0.1);
    else
        node = nextnode;
    end

end

end


function [dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff)
% pipeline consol

% failedflag = false;
if ~status.jobdone
    if echo_onoff, fprintf(' (%.2fsec)  failed in error code: %d\n', timecost, status.errorcode); end
    failedflag = true;
elseif status.errorcode==0
    if status.jobdone == 3
        if echo_onoff, fprintf(' (%.2fsec)  passed\n', timecost); end
    elseif status.jobdone == 5
        if echo_onoff, fprintf(' (%.2fsec)  skipped\n', timecost); end
    else
        if echo_onoff, fprintf(' (%.2fsec)  done\n', timecost); end
    end
    failedflag = false;
else
    if echo_onoff, fprintf(' (%.2fsec)  done but with error code: %d\n', timecost, status.errorcode); end
    failedflag = false;
end

if status.errorcode~=0
    if ischar(status.errormsg)
        warning(status.errormsg);
    elseif isa(status.errormsg, 'MException')
        warning(status.errormsg.getReport);
    else
        % error in throwing error?? what are you thinking?
        warning('Error in throwing error message!');
    end
    return;
end

[dataflow, status] = afternodedone(dataflow, status);

end


function [dataflow, status] = afternodedone(dataflow, status)
% run after every nodes in pipeline

if ~isfield(status.pipeline, status.nodename)
    % pass
    return
end

% status.jobdone=1:     normally done, to recycle and wake up next node, sleep (compatible with non-pineline nodes)
% status.jobdone=2:     technically done, to recycle and wake up next node, keep waking
% status.jobdone=3:     passed, sleep
% status.jobdone=4:     technically passed, to recycle, sleep
% status.jobdone=5:     skipping done, to recycle and wake up next node, sleep, just done but show up 'skipped'

% to wake up next node 
nextnode = status.pipeline.(status.nodename).nextnode;
switch status.jobdone
    case {1, 2, 5}
        % job done 
        if ~isempty(nextnode)
            % wake up next node
            status.pipeline.(nextnode).sleeping = false;
        end
    otherwise
        1;
end

% to recycle
if status.pipeline_onoff
    switch status.jobdone
        case {1, 2, 4, 5}
            % recycle
            if isfield(dataflow.pipepool, status.nodename)
                [dataflow.pipepool.(status.nodename), status.pipepool.(status.nodename)] = ...
                    pipepoolrecycle(dataflow.pipepool.(status.nodename), status.pipepool.(status.nodename));
            end
        otherwise
            1;
    end
end

% to sleep
switch status.jobdone
    case {1, 3, 4, 5}
        % go to sleep
        status.pipeline.(status.nodename).sleeping = true;
        % When we call a node functtion without pipeline updated, it will normally return status.jobdone=true. We shall let it
        % go to sleep, or the pipeline could be in endless looping.
    case 2
        % keep waking
        status.pipeline.(status.nodename).sleeping = false;
    otherwise
        % Note: the nodes can go to sleep in the node function.
        1;
end

% % is series done?
% prevnode = status.pipeline.(status.nodename).prevnode;
% if isempty(prevnode)
%     if status.pipeline.(status.nodename).sleeping
%         status.pipeline.(status.nodename).seriesdone = true;
%     end
% else
%     if status.pipeline.(status.nodename).sleeping && status.pipeline.(prevnode).seriesdone
%         status.pipeline.(status.nodename).seriesdone = true;
%     end
% end

% reset
status.currentjob = struct();

end


function [pooldata, poolstatus] = pipepoolrecycle(pooldata, poolstatus, recycle_onoff)

% poolfields
poolfields = fieldnames(pooldata);
if isempty(poolfields)
    return
end
% datalength
datalength = poolstatus.WritePoint - poolstatus.ReadPoint;

% recylestrategy
if nargin < 3
    switch poolstatus.recylestrategy
        case 0
            recycle_onoff = false;
        case 1
            % recycle
            recycle_onoff = true;
        case 2
            % recycle?
            recycle_onoff = poolstatus.WritePoint > poolstatus.warningstage;
        otherwise
            error('Pipeline consol error, unkown recylestrategy %d!', poolstatus.recylestrategy);
    end
end
% recycle
if recycle_onoff
    for ii = 1 : length(poolfields)
        if size(pooldata.(poolfields{ii}), 2) > datalength
            pooldata.(poolfields{ii})(:, 1:datalength) = ...
                pooldata.(poolfields{ii})(:, poolstatus.ReadPoint : poolstatus.WritePoint-1);
        elseif isstruct(pooldata.(poolfields{ii}))
            [pooldata.(poolfields{ii}), ~] = pipepoolrecycle(pooldata.(poolfields{ii}), poolstatus, recycle_onoff);
        end
    end
    poolstatus.ReadPoint = 1;
    poolstatus.WritePoint = datalength + 1;
end

end


function nodename = rollpoling(pipeline)

nodes = fieldnames(pipeline);
nodename = '';
priority0 = -99;
for ii = 1:length(nodes)
    if ~pipeline.(nodes{ii}).sleeping  && pipeline.(nodes{ii}).priority > priority0
        nodename = nodes{ii};
        priority0 = pipeline.(nodes{ii}).priority;
    end
end

end