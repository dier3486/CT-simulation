function [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status)
% support node, backup data
% [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status);

% parameters set in pipe
backupprm = prmflow.pipe.(status.nodename);

% backup index
if isfield(backupprm, 'index')
    backup_index = backupprm.index;
else
    backup_index = 1;
end

% backup fields
if isfield(backupprm, 'field')
    backkfields = backupprm.field;
else
    backkfields = 'rawdata';
end

% backup cell
if ~isfield(dataflow, 'backup')
    dataflow.backup = cell(1, backup_index);
end

if ~iscell(backkfields)
    if isfield(dataflow, backkfields)
        dataflow.backup{backup_index}.(backkfields) = dataflow.(backkfields);
    end
else
    for ii = 1:length(backkfields)
        if isfield(dataflow, backkfields{ii})
            dataflow.backup{backup_index}.(backkfields{ii}) = dataflow.(backkfields{ii});
        end
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end