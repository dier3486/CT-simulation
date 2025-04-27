function [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status)
% prio-step of a pipeline node

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;
pipeprm = nodeprm.pipeline;

% ini the flag to run the poststep
status.currentjob.torunpoststep = true;

% to skip the priostep
if ~pipeprm.priostep
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;    
    status.currentjob.pipeline.jobdone = true;
    status.jobdone = true;
    return;
end

% prio step #1
% prio step of pipepool
[dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode), status.currentjob.pipeline] = ...
    poolpriostep(dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode), pipeprm, status.currentjob.pipeline);

% trace
if status.currentjob.pipeline.isshotstart && status.debug.pooltrace_onoff
    dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'priostep');
end

% prio step #1.5
% set status.jobdone by status.currentjob.pipeline.jobdone;
status.jobdone = status.currentjob.pipeline.jobdone;
% check error
if status.jobdone == 0
    status.errorcode = status.currentjob.pipeline.errorcode;
    status.errormsg = [sprintf('Piostep of node %s error due to ', nodename) status.currentjob.pipeline.errormsg];
    status.currentjob.topass = true;
    return;
elseif status.jobdone == 3 || status.jobdone == 6
    status.currentjob.topass = true;
end

% prio step #2
% next pool (redirect by carrynode)
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

% check the next pool in resize (while next node is not NULL)
if ~isempty(dataflow.pipepool.(carrynode))
    if isavail(dataflow.pipepool.(carrynode)(1).poolsize)
        poolsize = dataflow.pipepool.(carrynode)(1).poolsize;
        forcetoresize = true;
    else
        poolsize = max(status.currentjob.pipeline.Index_out);
        forcetoresize = false;
    end
    dataflow.pipepool.(carrynode)(1).data = poolresize(dataflow.pipepool.(carrynode)(1).data, poolsize, ...
        dataflow.pipepool.(carrynode)(1).datafields, dataflow.pipepool.(nodename)(1).data, forcetoresize);
end

% prio step #3
% copy data to next pool
if pipeprm.kernellevel == 0 && ~pipeprm.iscarried
    datatoread = status.currentjob.pipeline.readnumber;
    if ~isempty(dataflow.pipepool.(carrynode))
        [dataflow.pipepool.(carrynode)(1).data, writenum] = ...
            pooldatacopy(dataflow.pipepool.(nodename)(1), dataflow.pipepool.(nodename)(1).data, ...
            dataflow.pipepool.(nextnode)(1), dataflow.pipepool.(carrynode)(1).data, datatoread, [], true, pipeprm.outputmethod);
        if datatoread ~= writenum
            status.jobdone = 0;
            status.errorcode = -301;
            status.errormsg = 'The data copy in prio-step was stucked!';
            status.currentjob.topass = true;
        end
    end
end

end
