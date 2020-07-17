function [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status)
% virus node, backup data
% [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status);

% parameters set in pipe
backupprm = prmflow.pipe.(status.nodename);

% backup index
if isfield(backupprm, 'index') && ~isempty(backupprm.index) 
    bkindex = num2str(backupprm.index);
else
    bkindex = num2str(status.seriesindex);
end

% dataflow, prmflow and status
if isfield(backupprm, 'dataflow')
    dataflow = backupdata(dataflow, backupprm.dataflow, bkindex);
end
if isfield(backupprm, 'prmflow')
    prmflow = backupdata(prmflow, backupprm.prmflow, bkindex);
end
if isfield(backupprm, 'status')
    status = backupdata(status, backupprm.status, bkindex);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function data = backupdata(data, bkfields, index)
% copy data.bkfields to data.bkfields_bk

% ext
name_ext = ['_bk' index];
% bkfields string to cell
if ~iscell(bkfields)
    bkfields = regexp(regexprep(bkfields, ' ', ''), ',', 'split');
end

% loop the fields to back up
for ii = 1:length(bkfields)
    if isfield(data, bkfields{ii})
        data.([char(bkfields{ii}) name_ext]) = data.(bkfields{ii});
    end
end
    
end
