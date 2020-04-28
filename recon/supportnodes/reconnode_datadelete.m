function [dataflow, prmflow, status] = reconnode_datadelete(dataflow, prmflow, status)
% virus node, delete data
% [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status);

% parameters set in pipe
backupprm = prmflow.pipe.(status.nodename);

% dataflow, prmflow and status
if isfield(backupprm, 'dataflow')
    dataflow = deletedata(dataflow, backupprm.dataflow);
end
if isfield(backupprm, 'prmflow')
    prmflow = deletedata(prmflow, backupprm.prmflow);
end
if isfield(backupprm, 'status')
    status = deletedata(status, backupprm.status);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function data = deletedata(data, delfields)
% copy data.bkfields to data.bkfields_bk

% delfields string to cell
if ~iscell(delfields)
    delfields = regexp(regexprep(delfields, ' ', ''), ',', 'split');
end
% isfield?
delfields = delfields(isfield(data, delfields));
% delete the fields
data = rmfield(data, delfields);
    
end
