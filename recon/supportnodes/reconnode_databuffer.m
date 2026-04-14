function [dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status)
% support node, to buffer data
% [dataflow, prmflow, status] = reconnode_databuffer(dataflow, prmflow, status);

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_databufferprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;
if ~pipeline_onoff
    % it shall be closed in prepare.
    error('reconnode_databuffer can only run in pipeline mode!');
end

% all the previous node sleeping?
nodesindex = 0 : status.pipeline.(nodename).index - 1;
allsleeping = nodestatecheck(status.pipeline, 'sleeping', nodesindex, 'all');

% prio-step
[dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
if status.currentjob.topass && ~allsleeping
    % error or pass
    return;
end

% main
% if to copy the buffer back to dataflow (after previous works done)
copytodataflow = nodeprm.copytodataflow;

% if to stuck the pipeline until previous works done, used to link a pipeline node with a non-pipeline node.
stuck_flag = nodeprm.stuck;

% copy data from currpool to buffer
Copynumber = status.currentjob.pipeline.readnumber; % it is not a temp variable
if dataflow.pipepool.(nodename)(1).iscarried
    inputcarrynode =  dataflow.pipepool.(nodename)(1).carrynode;
    inputcarryindex = dataflow.pipepool.(nodename)(1).carryindex;
else
    inputcarrynode = nodename;
    inputcarryindex = 1;
end
[dataflow.buffer.(nodename).pool(1).data, ~] = ...
    pooldatacopy(dataflow.pipepool.(nodename)(1), dataflow.pipepool.(inputcarrynode)(inputcarryindex).data, ...
            dataflow.buffer.(nodename).pool(1), dataflow.buffer.(nodename).pool(1).data, Copynumber, [], true);
% move the pointers in buffer pool
[~, dataflow.buffer.(nodename).pool(1)] = movepointsaftercopy([], dataflow.buffer.(nodename).pool(1), [], ...
    Copynumber, Copynumber);

% if to copy buffer data to dataflow
if copytodataflow && allsleeping
    % copy all the buffer data to dataflow (overwrite!)
    dataflow.buffer.(nodename).pool(1).ReadPoint = 1;
    writenum = dataflow.buffer.(nodename).pool(1).WritePoint - 1;
    targetpool = poolmirror(dataflow.buffer.(nodename).pool(1));
    targetpool.WritePoint = 1;
    [dataflow, ~] = pooldatacopy(dataflow.buffer.(nodename).pool(1), dataflow.buffer.(nodename).pool(1).data, ...
        targetpool, dataflow, writenum, [], true);

    % clear buffer
    dataflow.buffer.(nodename).pool(1).data = poolclear(dataflow.buffer.(nodename).pool(1).data);
    dataflow.buffer.(nodename).pool(1).ReadPoint = 1;
    dataflow.buffer.(nodename).pool(1).WritePoint = 1;
    dataflow.buffer.(nodename).pool(1).AvailPoint = 0;
    %     % series done
    %     status.pipeline.(nodename).seriesdone = true;
end

if stuck_flag
    % job done
    if allsleeping
        status.jobdone = 1;
        % 1: to wake up the next node which is supposed a non-pipeline node.
    else
        if Copynumber > 0
            status.jobdone = 4;
            % 4: technically pass, don't wake up the next node.
        else
            status.jobdone = 3;
            % 3: pass, did nothing
        end
    end
end

% post step
[dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);

end