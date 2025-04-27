function varargout = CTrecon(reconcfg, rawdatafile, offsetfile)
% CT reconstrcution
% images = CTrecon(reconxml);
% or
% [images, dataflow, prmflow] = CTrecon(reconxml, rawdatafile, offsetfile);

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

if ischar(reconcfg) || isstring(reconcfg)
    % try to read recon configure file (.xml or .json)
    reconcfg = readcfgfile(reconcfg);
end
if nargin<2
    rawdatafile = {};
end
if ~iscell(rawdatafile)
    rawdatafile = {rawdatafile};
end
Nraw = size(rawdatafile(:), 1);
if nargin<3
    offsetfile = {};
end
if ~iscell(offsetfile)
    offsetfile = {offsetfile};
end
Noffset = size(offsetfile(:), 1);

% series
if ~iscell(reconcfg.recon)
    reconcfg.recon = num2cell(reconcfg.recon);
end
Nseries = length(reconcfg.recon);

% clean path
if isfield(reconcfg, 'path')
    reconcfg = cleanpath(reconcfg, '', 'path');
    % copy the path to recon{ii}
    for ii = 1:Nseries
        reconcfg.recon{ii}.path = reconcfg.path;
    end
end

% ini outputs
images = cell(1, Nseries);
dataflow = struct();
prmflow = struct();

% ini status
status_ini = initialstatus(reconcfg.recon);

% loop the series
for iseries = 1:Nseries
    % initial status
    status = status_ini;
    status.seriesindex = iseries;

    % replace rawdata and offset
    if iseries<=Nraw
        status.reconcfg{iseries}.rawdata = char(rawdatafile{iseries});
    end
    if iseries<=Noffset
        status.reconcfg{iseries}.offset = char(offsetfile{iseries});
    end
    % recon access
    [dataflow, prmflow, status] = recon_access(dataflow, prmflow, status);
    
    % to return the images
    if isfield(dataflow, 'image')
        imagesize = prmflow.recon.imagesize;
        images{iseries} = reshape(dataflow.image, imagesize(2), imagesize(1), []);
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
