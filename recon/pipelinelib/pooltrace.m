function pipepool = pooltrace(pipepool, nodename, nextnode, opterator, branchnode, branchindex)
% recode the trace of the pools
%   dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, opterator);
% or,
%   dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, opterator, ...
%       {status.currentjob.pipeline.branchoutput(:).nextnode}, [status.currentjob.pipeline.branchoutput(:).poolindex]);

if ~isempty(nodename)
    for jj = 1:length(pipepool.(nodename))
        t = length(pipepool.(nodename)(jj).trace);
        pipepool.(nodename)(jj).trace(t+1).operator = [nodename ' ' opterator];
        pipepool.(nodename)(jj).trace(t+1) = poolmirror(pipepool.(nodename)(jj), ...
            pipepool.(nodename)(jj).trace(t+1));
    end
end

if ~isempty(nextnode)
    for jj = 1:length(pipepool.(nextnode))
        t = length(pipepool.(nextnode)(jj).trace);
        pipepool.(nextnode)(jj).trace(t+1).operator = [nodename ' ' opterator];
        pipepool.(nextnode)(jj).trace(t+1) = poolmirror(pipepool.(nextnode)(jj), ...
            pipepool.(nextnode)(jj).trace(t+1));
    end
end

if nargin >= 6 && ~isempty(branchnode)
    pipepool = branchpooltrace(pipepool, nodename, nextnode, opterator, branchnode, branchindex);
end

end

