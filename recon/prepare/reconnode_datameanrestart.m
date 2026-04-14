function [dataflow, prmflow, status] = reconnode_datameanrestart(dataflow, prmflow, status)
% it is an error 

status.errormsg = 'The datamean node can not be employed in restarting!';
status.jobdone = 0;

end