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

% go
echo_onoff = status.debug.echo_onoff;
if echo_onoff, fprintf('Recon Series %d\n', status.seriesindex); end

% initial steps (and initial GPU)
if echo_onoff, fprintf('  initial (GPU)...                           '); end
tic;
[dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');
timecost = toc;
[dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
if failedflag
    return;
end

% load calibration tables
if echo_onoff, fprintf('  load calibration tables...                 '); end
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
while ~status.seriesdone
    % wakeup the current node (it was returned by function pipelineconsol)
    node = status.nodename;
    % echo
    if echo_onoff, fprintf('  [recon node] %-30s', [node, ' ...']); end
    tic;

    % node main
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, node);
    timecost = toc;
    1;
    [dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff);
    if failedflag
        return;
    end
    
    % pipeline restart
    if status.torestart
        orignode = status.nodename;
        if echo_onoff, fprintf('  pipeline restarting...'); end
        [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'pipelinerestart');
        status.nodename = orignode;
    end
    
    % next node
    if strcmp(node, status.nodename)
        pause(0.1);
    end

end
% pipeline completed without error

end


function [dataflow, status, failedflag] = pipelineconsol(dataflow, status, timecost, echo_onoff)
% pipeline consol

% record the status.jobdone and walltime
if isfield(status.pipeline, status.nodename)
    status.pipeline.(status.nodename).jobdone = status.jobdone;
    status.pipeline.(status.nodename).walltime = status.pipeline.(status.nodename).walltime + timecost;
end

% echo time cost
if echo_onoff, fprintf(' (%.2fsec)  ', timecost); end
% check and echo jobdone
if ~status.jobdone
    fprintf(2, 'failed in error code: %d\n', status.errorcode);
    if status.errorcode==0 && isempty(status.errormsg)
        if echo_onoff, fprintf(2, 'Did you forget to show jobdone after work?\n'); end
    end
    failedflag = true;
elseif status.errorcode==0
    switch status.jobdone
        case 1
            if echo_onoff, fprintf('done\n'); end
        case 2
            if echo_onoff, fprintf('partly done\n'); end
        case 3
            if echo_onoff, fprintf('passed\n'); end
        case 4
            if echo_onoff, fprintf('technically passed\n'); end
        case 5
            if echo_onoff, fprintf('skipped\n'); end
        case 6 
            if echo_onoff, fprintf(2, 'stucked\n'); end
        case 7
            if echo_onoff, fprintf('call me again\n'); end
        otherwise
            if echo_onoff, fprintf(2, 'done in unknown status %d\n', status.jobdone); end
    end
    failedflag = false;
else
    if echo_onoff, fprintf(2, 'done with error in: %d\n', status.errorcode); end
    failedflag = false;
end

% echo error(warning) message
if status.errorcode~=0 && echo_onoff
    if ischar(status.errormsg)
        if ~status.jobdone
            fprintf(2, 'Error: ');
        else
            fprintf(2, 'Warning: ');
        end
        fprintf(2, status.errormsg);
        fprintf('\n');
    elseif isa(status.errormsg, 'MException')
        fprintf(2, status.errormsg.getReport);
    else
        % error in throwing error?? what are you thinking?
        fprintf(2, 'Error in throwing error message!\n');
    end
end

% failed
if failedflag
    return;
end

% operators depending on status.jobdone
afternodedone();

% error in afternodedone
if ~status.jobdone
    if echo_onoff, fprintf(2, 'Pipelineconsol failed in error code: %d, %s\n', status.errorcode, status.errormsg); end
    failedflag = true;
    return;
end

% new next node
newnode = rollpoling(status.pipeline);
status.seriesdone = isempty(newnode) || strcmpi(newnode, 'NULL');

if ~status.seriesdone
    % reset status.currentjob for the coming node
    status.currentjob = struct();
    status.currentjob.pipeline_onoff = status.pipeline.(newnode).pipeline_onoff;
    if status.pipeline.(newnode).pipeline_onoff
        status.currentjob.pipeline = struct();
    end
    status.currentjob.topass = false;
    % reset nodename
    status.nodename = newnode;
    % reset status.jobdone and error message
    status.jobdone = 0;
    status.errorcode = 0;
    status.errormsg = '';
else
    status.nodename = 'NULL';
    % seriesdone
    seriesdone();
end
% done

    % function afternodedone()
    function jobdone = afternodedone()
        % run after every nodes in pipeline

        if ~isfield(status.pipeline, status.nodename)  || status.jobdone==0
            % pass
            return
        end

        % status.jobdone=1:     normally done, to recycle and wake up next node, sleep (compatible with non-pineline nodes)
        % status.jobdone=2:     partly done, to recycle and wake up next node, keep waking
        % status.jobdone=3:     passed, sleep
        % status.jobdone=4:     technically passed, to recycle, sleep
        % status.jobdone=5:     skipping done, to recycle and wake up next node, sleep, just done but show up 'skipped'
        % status.jobdone=6:     stucking (due to output pool), wake up next node, keep waking
        % status.jobdone=7:     technically done, but do not wake up next node, keep waking (will be called again)

        % to wake up next node
        nextnode = status.pipeline.(status.nodename).nextnode;
        switch status.jobdone
            case {1, 2, 5, 6}
                % job done
                if ~isempty(nextnode) && ~strcmpi(nextnode, 'NULL')
                    % wake up next node (But do not wake up the NULL node)
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
                        pipelinerecycle();
%                         dataflow.pipepool.(status.nodename) = poolrecycle(dataflow.pipepool.(status.nodename));
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
            case {2, 6, 7}
                % keep waking
                status.pipeline.(status.nodename).sleeping = false;
            otherwise
                % Note: the nodes can go to sleep in the node function.
                1;
        end

        % check stuck
        if status.jobdone == 3
            % in passing
            prevnode = status.pipeline.(status.nodename).prevnode;
            if strcmpi(prevnode, 'NULL') || status.pipeline.(prevnode).jobdone == 6
                % 3+6 means stucked
                status.jobdone = 0;
                status.errorcode = 109;
                status.errormsg = 'pipeline stucked!';
            end
        end

        % return
        jobdone = status.jobdone;
        % function afternodedone END
    end

    % function seriesdone()
    function seriesdone()
        % run while status.seriesdone==true
        if echo_onoff
            fprintf('Series Completed');
            nodes = fieldnames(status.pipeline);
            if status.pipeline_onoff && ~isempty(nodes)
                fprintf(', online time cost:\n');
                walltime = 0;
                for ii = 1:length(nodes)
                    if any(strcmpi(nodes{ii}, {'loadrawdata', 'readraw', 'pipelinestuck'}))
                        continue;
                    end
                    nodestr = [nodes{ii} ':'];
                    fprintf('  %-26s%6.2fsec\n', nodestr, status.pipeline.(nodes{ii}).walltime);
                    walltime = walltime + status.pipeline.(nodes{ii}).walltime;
                end
                fprintf('  Total:%+20s%6.2fsec\n', '', walltime);
            else
                fprintf('\n');
            end
        end
        % function seriesdone END
    end
    
    % function recycle
    function pipelinerecycle(nodename)
        % to recylce the current node and to recyle the previous nodes which are carried by the current node
        % trace_onoff
        trace_onoff = status.debug.pooltrace_onoff;
        if nargin < 1
            nodename = status.nodename;
        end
        % loop the pools 
        for ii = 1:length(dataflow.pipepool.(nodename))
            if ~dataflow.pipepool.(nodename)(ii).iscarried
                % recycle
                [dataflow.pipepool.(nodename)(ii), ~, removenumber] = poolrecycle(dataflow.pipepool.(nodename)(ii));
                % trace
                if trace_onoff
                    t = length(dataflow.pipepool.(nodename)(ii).trace);
                    dataflow.pipepool.(nodename)(ii).trace(t+1).operator = [nodename ' recycle'];
                    dataflow.pipepool.(nodename)(ii).trace(t+1) = poolmirror(dataflow.pipepool.(nodename)(ii), ...
                        dataflow.pipepool.(nodename)(ii).trace(t+1));
                end
                % recycle the carried pools
                carriages = dataflow.pipepool.(nodename)(ii).carriages;
                if ~isempty(carriages) && removenumber~=0
                    for carriednode = carriages
                        for jj = 1:length(dataflow.pipepool.(carriednode{1}))
    % to find out the carried pool
    if strcmp(dataflow.pipepool.(carriednode{1})(jj).carrynode, nodename) && ...
            dataflow.pipepool.(carriednode{1})(jj).carryindex == ii
        % recycle
        dataflow.pipepool.(carriednode{1})(jj) = ...
            poolrecycle(dataflow.pipepool.(carriednode{1})(jj), struct(), removenumber);
        % trace
        if trace_onoff
            t = length(dataflow.pipepool.(carriednode{1})(jj).trace);
            dataflow.pipepool.(carriednode{1})(jj).trace(t+1).operator = [nodename ' recycle'];
            dataflow.pipepool.(carriednode{1})(jj).trace(t+1) = poolmirror(dataflow.pipepool.(carriednode{1})(jj), ...
                dataflow.pipepool.(carriednode{1})(jj).trace(t+1));
        end
    end
                        end
                    end
                end
            end
        end
        
%         % previous node
%         prevnode = status.pipeline.(nodename).prevnode;
%         if strcmpi(prevnode, 'NULL') || isempty(prevnode)
%             return;
%         else
%             % to recylce the previous node by a recursing,
%             % while the previous node's carrynode is the status.nodename it will be recylced and recursed again.
%             pipelinerecycle(prevnode);
%         end
    end
    % function recycle END
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