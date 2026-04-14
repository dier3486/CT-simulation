function [dataflow, prmflow, status] = helicalBPpriostep_dump(dataflow, prmflow, status)
% Helical BP priostep

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;
pipeprm = nodeprm.pipeline;

imagenumber = double(prmflow.recon.imagenumber);
% Nview = double(prmflow.recon.Nview);
Nviewskip = double(prmflow.recon.Nviewskip);
Nviewact = double(prmflow.recon.Nviewact);
viewbyimages = prmflow.recon.viewbyimages;

% iteration
if isfield(prmflow.recon, 'iteration_onoff')
    iteration_onoff = prmflow.recon.iteration_onoff;
else
    iteration_onoff = false;
end

% ini the flag to run the poststep
status.currentjob.torunpoststep = true;

% prio step #1, the prio steps of the pipepools
% pass due to nextpool is stucked
if dataflow.pipepool.(nextnode)(1).WriteStuck
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 6;
    % done and keep waking
    return;
end

% prio step #2, the next node and carry node
status.currentjob.nextnode = nextnode;
if ~isempty(dataflow.pipepool.(nextnode))
    if dataflow.pipepool.(nextnode)(1).iscarried
        carrynode = dataflow.pipepool.(nextnode)(1).carrynode;
    else
        carrynode = nextnode;
    end
else
    carrynode = nextnode;
end
status.currentjob.carrynode = carrynode;

% I know the nodetype is H-H.1.G

% check shot start 1
status.currentjob.pipeline.isshotstart = dataflow.pipepool.(nextnode)(1).isshotstart;
% close it
% dataflow.pipepool.(nextnode)(1).isshotstart = false;
% Isshotstart1 = dataflow.pipepool.(nodename)(1).ReadPoint == dataflow.pipepool.(nodename)(1).ReadStart;

if status.currentjob.pipeline.isshotstart
    % reset currpool.recycled
    dataflow.pipepool.(nodename)(1).recycled = false;
    % ini next pool
    dataflow.pipepool.(nextnode)(1).ReadStart = 1;
    dataflow.pipepool.(nextnode)(1).ReadEnd = imagenumber;
    dataflow.pipepool.(nextnode)(1).WriteStart = 1;
    dataflow.pipepool.(nextnode)(1).WriteEnd = imagenumber;
    dataflow.pipepool.(nextnode)(1).ReadPoint = 1;
    dataflow.pipepool.(nextnode)(1).WritePoint = 1;
%     % close the isshotstart
%     dataflow.pipepool.(nextnode)(1).isshotstart = false;
    % view skip
    dataflow.pipepool.(nodename)(1).ReadPoint = dataflow.pipepool.(nodename)(1).ReadStart + Nviewskip;
    dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd - Nviewskip;
    % Note: we moved the ReadPoint and ReadEnd but not the ReadStart
    % branches
    if iteration_onoff
        viewsubspace = prmflow.iteration.viewsubspace;
        subviewshift = prmflow.iteration.subviewshift;
        for ii = 1: viewsubspace
            % the branch node
            branchnode = status.pipeline.(nodename).branchnextnodes{ii};
            branchpoolindex = status.pipeline.(nodename).branchnextpoolindex(ii);
            if dataflow.pipepool.(branchnode)(branchpoolindex).WriteStuck
                % stuck due to the branchpool
                status.currentjob.pipeline.readnumber = 0;
                status.currentjob.pipeline.writenumber = 0;
                status.currentjob.pipeline.newAvail = 0;
                status.jobdone = 6;
%                 % skip to close the isshotstart
%                 dataflow.pipepool.(nextnode)(1).isshotstart = true;
                % to pass the kernel function
                status.currentjob.topass = true;
                % withdraw the initial of next pool
                if status.currentjob.pipeline.isshotstart
%                     dataflow.pipepool.(nextnode)(1).isshotstart = true;
                    dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
                end
                return;
            end
            Nview_ii = double(floor((prmflow.recon.Nviewact - subviewshift(ii)) / viewsubspace) + 1);
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadEnd = Nview_ii;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteEnd = Nview_ii;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadPoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WritePoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).AvailPoint = 0;
        end
        1;
    end
end

% trace
if status.currentjob.pipeline.isshotstart && status.debug.pooltrace_onoff
    if ~iteration_onoff
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'priostep');
    else
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'priostep', ...
            status.pipeline.(nodename).branchnextnodes, status.pipeline.(nodename).branchnextpoolindex);
    end
end

% check shot start 2
% status.currentjob.pipeline.isshotstart = dataflow.pipepool.(nodename)(1).ReadPoint == ...
%     dataflow.pipepool.(nodename)(1).ReadStart + Nviewskip;

% check shot end 1 (the input data has reached the end)
Isshotend1 = dataflow.pipepool.(nodename)(1).WritePoint == dataflow.pipepool.(nodename)(1).WriteEnd + 1;

% minlimit/maxlimit
status.currentjob.pipeline.minlimit = pipeprm.inputminlimit;
status.currentjob.pipeline.maxlimit = pipeprm.inputmaxlimit;

% avail inputs
currAvailNumber = min(dataflow.pipepool.(nodename)(1).AvailPoint, dataflow.pipepool.(nodename)(1).ReadEnd) ...
    - dataflow.pipepool.(nodename)(1).ReadPoint + 1;
% I know, we have moved the ReadEnd therefore the AvailPoint could run over the ReadEnd.
n = min(pipeprm.inputmaxlimit, currAvailNumber);

% space left in nextpool
nextleft = dataflow.pipepool.(nextnode)(1).poolsize - dataflow.pipepool.(nextnode)(1).WritePoint + 1;

% m (active images)
% AvailV = [0  dataflow.pipepool.(nodename)(1).AvailPoint - dataflow.pipepool.(nodename)(1).ReadPoint] + ...
%     prmflow.recon.viewread + 1;
AvailV = [1 n] + prmflow.recon.viewread;
s = (viewbyimages(2, :) >= AvailV(1)) & (viewbyimages(1, :) <= AvailV(2));
m = sum(s);
if m > nextleft
    m = nextleft;
    n = viewbyimages(1, find(s, 1) + nextleft) - AvailV(1);
    if n < status.currentjob.pipeline.minlimit || nextleft <= 0
        % stucking
        status.currentjob.pipeline.readnumber = 0;
        status.currentjob.pipeline.writenumber = 0;
        status.currentjob.pipeline.newAvail = 0;
        status.jobdone = 6;
        % to pass the kernel function
        status.currentjob.topass = true;
        % withdraw the initial of next pool
        if status.currentjob.pipeline.isshotstart
%             dataflow.pipepool.(nextnode)(1).isshotstart = true;
            dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
        end
        return;
    end
elseif n < status.currentjob.pipeline.minlimit && ~Isshotend1
    % not enough input views
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 3;
    % to pass the kernel function
    status.currentjob.topass = true;
    % withdraw the initial of next pool
    if status.currentjob.pipeline.isshotstart
%         dataflow.pipepool.(nextnode)(1).isshotstart = true;
        dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
    end
    return;
end
if n < currAvailNumber
    % I know nextleft>0
    if currAvailNumber - n >= pipeprm.inputminlimit || Isshotend1
        % partly done
        status.jobdone = 2;
    else
        status.jobdone = 1;
    end
else
    % done
    status.jobdone = 1;
end

% check shot end 2 (will output the shot end)
status.currentjob.pipeline.isshotend = Isshotend1 && n == currAvailNumber;

% Index_in, Index_out
status.currentjob.pipeline.Index_in = ...
    [dataflow.pipepool.(nodename)(1).ReadPoint dataflow.pipepool.(nodename)(1).ReadPoint+n-1];
status.currentjob.pipeline.Index_out = ...
    [dataflow.pipepool.(nextnode)(1).WritePoint dataflow.pipepool.(nextnode)(1).WritePoint+m-1];
% Do2i
% status.currentjob.pipeline.Do2i = ...
%     dataflow.pipepool.(nodename)(1).ReadPoint - dataflow.pipepool.(nodename)(1).ReadStart - Nviewskip;
% % I know the Do2i will be used to calculate the Zview, not what it was.
status.currentjob.pipeline.Do2i = 0;
% Nexpand (has been count)
status.currentjob.pipeline.Nexpand = 0;

% readnumber
status.currentjob.pipeline.readnumber = n;

% writenumber
s = (viewbyimages(2, :) >= AvailV(1)) & (min(viewbyimages(2, :), Nviewact) <= AvailV(1)+n-1);
writenumber = sum(s);
status.currentjob.pipeline.writenumber = writenumber;
% newAvail
% status.currentjob.pipeline.newAvail = writenumber;
Navail = dataflow.pipepool.(nextnode)(1).WritePoint - dataflow.pipepool.(nextnode)(1).AvailPoint - 1 + writenumber;
if Navail < pipeprm.outputminlimit && ~status.currentjob.pipeline.isshotend
    status.currentjob.pipeline.newAvail = 0;
else
    status.currentjob.pipeline.newAvail = Navail;
end
% index_avail
status.currentjob.pipeline.Index_avail = ...
    dataflow.pipepool.(nextnode)(1).AvailPoint + [1 status.currentjob.pipeline.newAvail];
    
% fix jobdone
if status.currentjob.pipeline.newAvail == 0
    if status.jobdone == 1
        status.jobdone = 4;
    elseif status.jobdone == 2
        status.jobdone = 7;
    end
end

% branchs in status.currentjob
if iteration_onoff
    viewsubspace = prmflow.iteration.viewsubspace;
    status.currentjob.pipeline.branchoutput(viewsubspace) = struct();
    for isub = 1 : viewsubspace
        % the branch node
        branchnode = status.pipeline.(nodename).branchnextnodes{isub};
        branchpoolindex = status.pipeline.(nodename).branchnextpoolindex(isub);
        if dataflow.pipepool.(branchnode)(branchpoolindex).iscarried
            carrynode = dataflow.pipepool.(branchnode)(branchpoolindex).carrynode;
            carryindex = dataflow.pipepool.(branchnode)(branchpoolindex).carryindex;
        else
            carrynode = branchnode;
            carryindex = branchpoolindex;
        end
        status.currentjob.pipeline.branchoutput(isub).nextnode = branchnode;
        status.currentjob.pipeline.branchoutput(isub).poolindex = branchpoolindex;
        status.currentjob.pipeline.branchoutput(isub).carrynode = carrynode;
        status.currentjob.pipeline.branchoutput(isub).carryindex = carryindex;
        status.currentjob.pipeline.branchoutput(isub).writenumber = nan;    % to be set
        status.currentjob.pipeline.branchoutput(isub).newAvail = nan;       % to be set
        status.currentjob.pipeline.branchoutput(isub).Index = ...
            [dataflow.pipepool.(branchnode)(branchpoolindex).WritePoint -1];   % to be set, or skip
    end
end

end