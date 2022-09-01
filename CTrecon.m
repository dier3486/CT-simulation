function varargout = CTrecon(reconxml, rawdatafile)
% CT reconstrcution
% images = CTrecon(reconxml);
% or
% [images, dataflow, prmflow] = CTrecon(reconxml, rawdatafile);

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

if ischar(reconxml) || isstring(reconxml)
    % try to read recon xml file
    reconxml = readcfgfile(reconxml);
end
if nargin<2
    rawdatafile = {};
end
if ~iscell(rawdatafile)
    rawdatafile = {rawdatafile};
end
Nraw = size(rawdatafile(:), 1);

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
status.taskUID = dicomuid();
% loop the series
for iseries = 1:Nseries
    status.seriesindex = iseries;
    % replace rawdata
    if iseries<=Nraw
        status.reconcfg{iseries}.rawdata = char(rawdatafile{iseries});
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

% return images
varargout{1} = images;
% return pipe
varargout{2} = dataflow;
varargout{3} = prmflow;
varargout{4} = status;

end
