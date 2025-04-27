function scp = outputbranchpool(scp, pipepool, branchindex, currnode, poolindex, writenumber, newAvail, Index_out)
% the prio-step of the branch ouput pools, 
%   status.currentjob.pipeline = outputbranchpool(status.currentjob.pipeline, dataflow.pipepool, ...
%     branchindex, currnode, poolindex, writenumber, newAvail)
% the scp shall be the status.currentjob.pipeline

if isempty(poolindex)
    poolindex = branchindex + 1;
end
if nargin < 7 || isempty(newAvail)
    newAvail = writenumber;
end

if pipepool.(currnode)(poolindex).iscarried
    carrynode = pipepool.(currnode)(poolindex).carrynode;
    carryindex = pipepool.(currnode)(poolindex).carryindex;
else
    carrynode = currnode;
    carryindex = poolindex;
end

scp.branchoutput(branchindex).nextnode = currnode;
scp.branchoutput(branchindex).poolindex = poolindex;
scp.branchoutput(branchindex).carrynode = carrynode;
scp.branchoutput(branchindex).carryindex = carryindex;
scp.branchoutput(branchindex).writenumber = writenumber;
scp.branchoutput(branchindex).newAvail = newAvail;
if nargin < 8 || isempty(Index_out)
    scp.branchoutput(branchindex).Index = ...
        [pipepool.(currnode)(poolindex).WritePoint  pipepool.(currnode)(poolindex).WritePoint + writenumber - 1];
else
    Index_out = Index_out - Index_out(1) + pipepool.(currnode)(poolindex).WritePoint;
    scp.branchoutput(branchindex).Index = Index_out;
end

end