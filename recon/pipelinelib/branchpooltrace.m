function pipepool = branchpooltrace(pipepool, nodename, nextnode, opterator, branchnode, branchindex)
% recode the trace of the pools
%   dataflow.pipepool = branchpooltrace(pipepool, nodename, nextnode, opterator, branchnode, branchindex)


for ii = 1 : length(branchnode)
    if isempty(branchnode{ii})
        continue;
    end
    if ~strcmp(branchnode{ii}, nextnode)
        t = length(pipepool.(branchnode{ii})(branchindex(ii)).trace);
        pipepool.(branchnode{ii})(branchindex(ii)).trace(t+1).operator = [nodename ' ' opterator];
        pipepool.(branchnode{ii})(branchindex(ii)).trace(t+1) = poolmirror(pipepool.(branchnode{ii})(branchindex(ii)), ...
            pipepool.(branchnode{ii})(branchindex(ii)).trace(t+1));
    end
end


