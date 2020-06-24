function varargout = CRISrecon(reconxml, rawdatafile)
% CT reconstrcution
% images = CTrecon(reconxml);
% or
% [images, dataflow, prmflow] = CTrecon(reconxml, rawdatafile);

% clean inputs
% rawdatafile
if nargin<2
    rawdatafile = {};
end
if ~iscell(rawdatafile)
    rawdatafile = {rawdatafile};
end
Nraw = size(rawdatafile(:), 1);
% reconxml
if isempty(reconxml)
    reconxml = struct();
    reconxml.recon = cell(1, Nraw);
    reconxml.recon(:) = struct();
elseif ischar(reconxml) || isstring(reconxml)
    % try to read recon xml file
    reconxml = readcfgfile(char(reconxml));
end

% series
if ~iscell(reconxml.recon)
    reconxml.recon = {reconxml.recon};
end
Nseries = length(reconxml.recon);

% ini outputs
images = cell(1, Nseries);
dataflow = struct();
prmflow = struct();

% ini status
status = struct();
status.reconcfg = reconxml.recon;
% loop the series
for iseries = 1:Nseries
    status.seriesindex = iseries;
    % replace rawdata
    if iseries<=Nraw
        status.reconcfg{iseries}.rawdata = char(rawdatafile{iseries});
    end
    status.reconcfg{iseries}.rawdata = checkrawfile(status.reconcfg{iseries}.rawdata, '.(raw|pd)');
    % CRIS supports
    rawxml = CRISgetrawxml(status.reconcfg{iseries}.rawdata);
    if ~isempty(rawxml)
        prmflow.external.rawxml = rawxml;
    end
    
    % recon access
    [dataflow, prmflow, status] = recon_access(status, 1, dataflow, prmflow);
    
    % to return the images
    if isfield(dataflow, 'image')
        images{iseries} = dataflow.image;
    end
    % NOTE: if we set to output the image to dicom, e.g. recon.pipe.dataoutput.files = 'dicomimage_namekey', the saved file's
    % name will be returned in prmflow.output.dicomimage
end

% return
varargout{1} = images;
varargout{2} = dataflow;
varargout{3} = prmflow;

end


