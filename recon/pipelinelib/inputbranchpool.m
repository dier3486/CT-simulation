function scp = inputbranchpool(scp, pipepool, branchindex, currnode, poolindex, readnumber)
% the prio-step of the branch input pools, 
%   status.currentjob.pipeline = inputbranchpool(status.currentjob.pipeline, dataflow.pipepool, branchindex, currnode, [], ...
%     readnumber);
% the scp shall be the status.currentjob.pipeline

if isempty(poolindex)
    poolindex = branchindex + 1;
end
if pipepool.(currnode)(poolindex).iscarried
    carrynode = pipepool.(currnode)(poolindex).carrynode;
    carryindex = pipepool.(currnode)(poolindex).carryindex;
else
    carrynode = currnode;
    carryindex = poolindex;
end

scp.branchinput(branchindex).poolindex = poolindex;
scp.branchinput(branchindex).carrynode = carrynode;
scp.branchinput(branchindex).carryindex = carryindex;
scp.branchinput(branchindex).readnumber = readnumber;
scp.branchinput(branchindex).Index = ...
    [pipepool.(currnode)(poolindex).ReadPoint  pipepool.(currnode)(poolindex).ReadPoint + readnumber - 1];
% scp.readnumber(poolindex) = readnumber;

end