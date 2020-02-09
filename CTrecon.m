function varargout = CTrecon(reconxml)
% CT reconstrcution
% output = CTrecon(reconxml);

if ischar(reconxml)
    % try to read recon xml file
    reconxml = readcfgfile(reconxml);
end

% series
if ~iscell(reconxml.recon)
    reconxml.recon = {reconxml.recon};
end
Nseries = length(reconxml.recon);

% ini outputs
images = cell(1, Nseries);
if nargout>1
    dataflow = struct();
    prmflow = struct();
end

status = struct();
status.reconcfg = reconxml.recon;
% loop the series
for iseries = 1:Nseries
    status.series_index = iseries;
    % recon access
    if nargout>1
        [dataflow, prmflow, status] = recon_access(status, 1, dataflow, prmflow);
    else
        [dataflow, prmflow, status] = recon_access(status, 1);
    end
    
    % to return the images
    if isfield(dataflow, 'image')
        images{iseries} = dataflow.image;
    end
end

% return
varargout{1} = images;
varargout{2} = dataflow;
varargout{3} = prmflow;

end
