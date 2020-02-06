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
Nseries = length({reconxml.recon});

% ini outputs
images = cell(1, Nseries);

status = struct();
% loop the series
for iseries = 1:length(Nseries)
    status.reconcfg = reconxml.recon{iseries};
    status.series_index = iseries;
    % recon access
    [dataflow, prmflow, status] = recon_access(status, 1);
    
    % to return the images
    if isfield(dataflow, 'image')
        images{iseries} = dataflow.image;
    end
end

% return
varargout{1} = images;
varargout{2} = prmflow;

end
