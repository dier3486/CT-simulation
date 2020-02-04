function varargout = recon_access(reconxml)
% recon & cali governing function
% output = recon_access(reconxml);

if ischar(reconxml)
    % try to read recon xml file
    reconxml = readcfgfile(reconxml);
end

% series
if ~iscell(reconxml.recon)
    reconxml.recon = {reconxml.recon};
end
Nseries = length({reconxml.recon});

% ini outputs
images = cell(1, Nseries);

% loop the series
for iseries = 1:length(Nseries)
    status.reconcfg = reconxml.recon{iseries};
    status.series_index = iseries;
    % initial
    dataflow = struct();
    prmflow = struct();
    % initial steps
    fprintf('Recon Series %d\n', iseries);
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'initial');
    fprintf('  load calibration tables...');
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadcorrs');
    fprintf(' done\n');
    
    % load rawdata
    fprintf('  read rawdata...');
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, 'loadrawdata');
    fprintf(' done\n');
    % for large data we should employ view buffer in loading data, TBC
    
    % run pipe nodes
    pipefields = fieldnames(prmflow.pipe);
    for i_node = 1:length(pipefields)
        node = pipefields{i_node};
        fprintf('  [recon node] %s...', node);
        [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, node);
        fprintf(' done\n');
    end
    fprintf('Done\n');
    
    % to return the images
    if isfield(dataflow, 'image')
        images{iseries} = dataflow.image;
    end
end

% return
varargout{1} = images;

end
