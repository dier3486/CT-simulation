function [dataflow, prmflow, status] = reconnode_sleep(dataflow, prmflow, status)
% support node, sleep
% [dataflow, prmflow, status] = reconnode_sleep(dataflow, prmflow, status);

% parameters set in pipe
sleepprm = prmflow.pipe.(status.nodename);

if isfield(sleepprm, 'time')
    pause();
else
    pause(sleepprm.time);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end