function [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status)
% post-step of a pipeline node

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;

if isfield(status.currentjob, 'torunpoststep')
    torunpoststep = status.currentjob.torunpoststep;
else
    torunpoststep = true;
end

% to skip the poststep
if ~nodeprm.pipeline.poststep || ~torunpoststep
    return;
end

% numbers
readnumber = double(status.currentjob.pipeline.readnumber);
writenumber = double(status.currentjob.pipeline.writenumber);
newAvail = double(status.currentjob.pipeline.newAvail);
Nexpand = double(status.currentjob.pipeline.Nexpand);

% move the points: currpool.ReadPoint, nextpool.WritePoint and nextpool.AvailPoint
[dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode)] = ...
    movepointsaftercopy(dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode), readnumber, writenumber, newAvail);
% the readnumber, writenumber and newAvail can be vector to move the branch pools
% or we will check the status.currentjob.pipeline.branchinput and .branchoutput to do that.

% check overflow
if ~dataflow.pipepool.(nodename)(1).circulatemode
    if dataflow.pipepool.(nodename)(1).ReadPoint > dataflow.pipepool.(nodename)(1).poolsize + 1
        status.errorcode = 310;
        status.errormsg = sprintf('The pool %s(1) overflowed in reading.', nodename);
    end
end
if ~isempty(dataflow.pipepool.(nextnode)) && ~dataflow.pipepool.(nextnode)(1).circulatemode
    if dataflow.pipepool.(nextnode)(1).WritePoint + Nexpand > dataflow.pipepool.(nextnode)(1).poolsize + 1
        status.errorcode = 310;
        status.errormsg = sprintf('The pool %s(1) overflowed in writing.', nextnode);
    end
end

% branch input pools
if isfield(status.currentjob.pipeline, 'branchinput') && ~isempty(fieldnames(status.currentjob.pipeline.branchinput))
    for ibranch = 1 : length(status.currentjob.pipeline.branchinput)
        branchindex = status.currentjob.pipeline.branchinput(ibranch).poolindex;
        if isempty(branchindex)
            continue;
        end
        if branchindex <= length(readnumber) && readnumber(branchindex) > 0
            % the readnumber has been set
            continue;
        end
        branchreadnumber = double(status.currentjob.pipeline.branchinput(ibranch).readnumber);
        dataflow.pipepool.(nodename)(branchindex) = ...
            movepointsaftercopy(dataflow.pipepool.(nodename)(branchindex), [], branchreadnumber);
        % check overflow
        if ~dataflow.pipepool.(nodename)(branchindex).circulatemode
            if dataflow.pipepool.(nodename)(branchindex).ReadPoint > dataflow.pipepool.(nodename)(branchindex).poolsize + 1
                status.errorcode = 311;
                status.errormsg = sprintf('The pool %s(%d) overflowed in reading.', nodename, branchindex);
            end
        end
    end
end

% branch output pools
if isfield(status.currentjob.pipeline, 'branchoutput') && ~isempty(fieldnames(status.currentjob.pipeline.branchoutput))
    for ibranch = 1 : length(status.currentjob.pipeline.branchoutput)
        branchnode = status.currentjob.pipeline.branchoutput(ibranch).nextnode;
        branchindex = status.currentjob.pipeline.branchoutput(ibranch).poolindex;
        if isempty(branchnode)
            continue;
        end
        if strcmp(branchnode, nextnode)
            if branchindex <= length(writenumber) && writenumber(branchindex) > 0
                % the branchnode is just the nextnode and the writenumber has been set
                continue;
            elseif branchindex <= length(newAvail) && newAvail(branchindex) > 0
                % the branchnode is just the nextnode and newAvail has been set
                continue;
            end
        end % else the branchnode is not the nextnode
        branchwritenumber = double(status.currentjob.pipeline.branchoutput(ibranch).writenumber);
        brnachnewAvial = double(status.currentjob.pipeline.branchoutput(ibranch).newAvail);
        [~, dataflow.pipepool.(branchnode)(branchindex)] = ...
            movepointsaftercopy([], dataflow.pipepool.(branchnode)(branchindex), [], branchwritenumber, brnachnewAvial);
        % check overflow
        if ~dataflow.pipepool.(branchnode)(branchindex).circulatemode
            if dataflow.pipepool.(branchnode)(branchindex).WritePoint > dataflow.pipepool.(branchnode)(branchindex).poolsize + 1
                status.errorcode = 311;
                status.errormsg = sprintf('The pool %s(%d) overflowed in writing.', branchnode, branchindex);
            end
        end
    end
end

% trace
if status.debug.pooltrace_onoff
    if isfield(status.currentjob.pipeline, 'branchoutput') && ~isempty(fieldnames(status.currentjob.pipeline.branchoutput))
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'poststep', ...
            {status.currentjob.pipeline.branchoutput(:).nextnode}, [status.currentjob.pipeline.branchoutput(:).poolindex]);
    else
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'poststep');
    end
end

end