function [dataflow, prmflow, status] = reconinitial(dataflow, prmflow, status)
% recon initial
% [dataflow, prmflow, status] = reconinitial(dataflow, prmflow, status)

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

% copy status.reconcfg to prmflow
reconcfg = status.reconcfg{status.seriesindex};

if isempty(reconcfg)
    % empty configure?
    status.jobdone = false;
    status.errorcode = 1;
    status.errormsg = '[reconinitial] empty recon configure';
    return;
end

% copy reconcfg to prmflow
prmflow = structmerge(reconcfg, prmflow, 0, 0);
% but maintain the extra fields in prmflow

% reload sub-config file
prmflow = subconfigure(prmflow);

% external supports
if isfield(prmflow, 'external')
    switch prmflow.external.Manufacturer
        case 'SINOVISION'
            % CRIS platform of SINOVISION
            prmflow = CRIS2prmflow(prmflow, prmflow.external.rawxml);
        otherwise
            1;
    end
end

% clean prmflow
prmflow = iniprmclean(prmflow);

% ini GPU
if isfield(prmflow.system, 'GPUdeviceindex')
    status.GPUinfo = initialGPU(prmflow.system.GPUdeviceindex);
    % I know when the GPUdeviceindex=0 the GPUinfo is []
else
    status.GPUinfo = [];
end

% series UID
status.seriesUID = dicomuid();

% datablock and pipeline initial
[dataflow, prmflow, status] = pipelineinitial(dataflow, prmflow, status);

% while (status.seriesdone) the pipeline will be closed
status.seriesdone = false;

% while (status.torestart) the pipeline will hold a restarting process
status.torestart = false;

% nodes can use the status.currentjob to save inner parameters
status.currentjob = struct();

% to return the status of current node
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function prmflow = iniprmclean(prmflow)
% to fill up the paramters which could be used in recon but not configured
% hard code

% clean path
if isfield(prmflow, 'path')
    prmflow = cleanpath(prmflow, '', 'path');
end

% in system
if isfield(prmflow, 'system') && ~isempty(prmflow.system)
    % collimatorexplain
    if ~isfield(prmflow.system, 'collimatorexplain') || isempty(prmflow.system.collimatorexplain)
        prmflow.system.collimatorexplain = [];
    elseif ischar(prmflow.system.collimatorexplain)
        prmflow.system.collimatorexplain = readcfgfile(prmflow.system.collimatorexplain);
        if ~iscell(prmflow.system.collimatorexplain.collimator)
            prmflow.system.collimatorexplain.collimator = num2cell(prmflow.system.collimatorexplain.collimator);
        end
    else
        error('Illegal prmflow.system.collimatorexplain class: %s!', class(prmflow.system.collimatorexplain));
    end
    % the collimatorexplain is used to exlpaint the protocol' collimator settings which are strings, e.g. '16x0.5'.
    % filetagsrule
    if ~isfield(prmflow.system, 'filetagsrule')
        prmflow.system.filetagsrule = struct();
    elseif ischar(prmflow.system.filetagsrule)
        prmflow.system.filetagsrule = readcfgfile(prmflow.system.filetagsrule);
    else
        error('Illegal prmflow.system.filetagsrule class: %s!', class(prmflow.system.filetagsrule));
    end
    % the filetagsrule is used to sparse or produce the filenames by the protocol, e.g. to present the collimator setting in
    % filename.
    % filematchrule
    if ~isfield(prmflow.system, 'filematchrule')
        prmflow.system.filematchrule = struct();
    elseif ischar(prmflow.system.filematchrule)
        prmflow.system.filematchrule = readcfgfile(prmflow.system.filematchrule);
    else
        error('Illegal prmflow.system.filematchrule class: %s!', class(prmflow.system.filematchrule));
    end
    % the filematchrule is used to match the fils (mostly calibration tables) with the protocol.
    % GPU+
    if ~isfield(prmflow.system, 'GPUdeviceindex')
        if gpuDeviceCount > 0
            prmflow.system.GPUdeviceindex = 1;
        else
            prmflow.system.GPUdeviceindex = 0;
        end
    end
    % GPUdeviceindex is (are) the device numbers of the GPU(s).
    % set dicom dictionary
    if isfield(prmflow.system, 'dicomdictionary') && ~isempty(prmflow.system.dicomdictionary)
        dicomdict('set', prmflow.system.dicomdictionary);
    end
    % the dicomdictionary is used to sparse or package the private DICOM files.
    % maxFOV
    if ~isfield(prmflow.system, 'maxFOV')
        % default maxFOV
        prmflow.system.maxFOV = 500.1;
        warning('Missing the maxFOV in recon system configure!');
    end
    % max scan/recon FOV, have to have.
    % default poolsize
    if ~isfield(prmflow.system, 'defaultrawpoolsize')
        % default raw pool size
        prmflow.system.defaultrawpoolsize = 512;
    end
    % defaultrawpoolsize is the default poolsize of rawdata used in pipeline-on mode
    if ~isfield(prmflow.system, 'defaultimagepoolsize')
        % default image pool size
        prmflow.system.defaultimagepoolsize = 128;
    end
    % defaultimagepoolsize is the default poolsize of images used in pipeline-on mode
    % max view number
    if ~isfield(prmflow.system, 'maxviewnumber')
        % default maxviewnumber
        prmflow.system.maxviewnumber = inf;
    end
    % maxviewnumber is the threshold of the view number restarting. (But NOT the max limit of the viewnumber, whose limit shall
    % be infinite.)
    % reshape of focalposition
    if isfield(prmflow.system, 'focalposition')
        prmflow.system.focalposition = reshape(prmflow.system.focalposition, [], 3);
    end
else
    prmflow.system = struct();
end

% protocol
% rawdata path/filename
if isfield(prmflow, 'rawdata')
    % copy rawdata path to .protocol
    prmflow.protocol.rawdata = prmflow.rawdata;     
elseif ~isfield(prmflow.protocol, 'rawdata')
    prmflow.protocol.rawdata = '';
end
% output path
if isfield(prmflow, 'outputpath')
    % copy rawdata outputpath to .protocol
    prmflow.protocol.outputpath = prmflow.outputpath;
end
if ~isfield(prmflow.protocol, 'outputpath') || isempty(prmflow.protocol.outputpath)
    % default outputpath
    filepath = fileparts(prmflow.protocol.rawdata);
    if isempty(filepath)
        % default outputpath is .\
        prmflow.protocol.outputpath = '.';
    else
        % default outputpath is the path of rawdata + \dcm
        prmflow.protocol.outputpath = fullfile(filepath, 'dcm');
    end
end
% empty prmflow.pipe
if ~isfield(prmflow, 'pipe') || isempty(prmflow.pipe)
    prmflow.pipe = struct();
end
prmflow.protocol = iniprotocolclean(prmflow.protocol, prmflow.pipe, prmflow.system);

% pipeline replicate
if isempty(prmflow.protocol.pipelinereplicate)
    % default pipelinereplicate depends on the datablock
    if isempty(prmflow.protocol.datablock)
        prmflow.protocol.pipelinereplicate = false;
    else
        prmflow.protocol.pipelinereplicate = true;
    end
end

% datablock
if prmflow.protocol.pipelinereplicate && isempty(prmflow.protocol.datablock) && isfield(prmflow.system, 'rawdatablock')
    % pipeline_on but missed to set datablock
    prmflow.protocol.datablock = prmflow.system.rawdatablock;
    % I know the protocol.datablock has been set to [] if it wasn't set in recon xml.
end


% IOstandard
if ~isfield(prmflow, 'IOstandard')
    prmflow.IOstandard = [];
end
% ini corrtable (always)
% if ~isfield(prmflow, 'corrtable')
%     prmflow.corrtable = struct();
% end
prmflow.corrtable = struct();
% ini raw, correction, rebin, recon and image in prmflow (always)
prmflow.raw = struct();
prmflow.correction = struct();
prmflow.rebin = struct();
prmflow.recon = struct();
prmflow.image = struct();
% NOTE: sometimes we need to maintain data in prmflow for follow-up series, so we don't clean all of the informations in 
% prmflow, so carefully in using that will be cleaned fields.

% pipe
if ~isfield(prmflow, 'pipe')
    prmflow.pipe = struct();
end
for pipenodes = fieldnames(prmflow.pipe)'
    if iscell(prmflow.pipe.(pipenodes{1}))
        warning('The pipeline nodes'' name is nonreusable!');
        prmflow.pipe.(pipenodes{1}) = prmflow.pipe.(pipenodes{1}){1};
    end
    if ~isstruct(prmflow.pipe.(pipenodes{1}))
        % replace the empty settings to struct();
        prmflow.pipe.(pipenodes{1}) = struct();
    end
end

end


function protocol = iniprotocolclean(protocol, pipe, system)
% to fill up the paramters in protocol
% hard code

% imagethickness
if isfield(protocol, 'collimator')
    collitoken = regexp(protocol.collimator, 'x([\d \.]+)\>', 'tokens');
else
    collitoken = [];
end
if isempty(collitoken)
    % what??
    CollimatedSliceThickness = 0;
else
    CollimatedSliceThickness = str2double(collitoken{1}{1});
end
if ~isfield(protocol, 'imagesize')
    % default imagesize 512x512
    protocol.imagesize = [512, 512];
elseif length(protocol.imagesize(:)) == 1
    protocol.imagesize = [protocol.imagesize, protocol.imagesize];
else
    protocol.imagesize = protocol.imagesize(:)';
end
if ~isfield(protocol, 'imagethickness')
    protocol.imagethickness = CollimatedSliceThickness;
    % set imagethickness with CollimatedSliceThickness
end
% imageincrement
if ~isfield(protocol, 'imageincrement')
    protocol.imageincrement = protocol.imagethickness;
end
% couchdirection
if ~isfield(protocol, 'couchdirection')
    if isfield(protocol, 'shotcouchstep')  && isfield(protocol, 'couchspeed')
        protocol.couchdirection = sign(sign(protocol.shotcouchstep) + sign(protocol.couchspeed));
    else
        protocol.couchdirection = 0;
    end
end
% gantrytilt
if ~isfield(protocol, 'gantrytilt')
    protocol.gantrytilt = 0;
    % in 180 degree
end
% gantrytilt
if ~isfield(protocol, 'focalspot')
    protocol.focalspot = 0;
end
% reconcenter
if ~isfield(protocol, 'reconcenter')
    protocol.reconcenter = [0 0];
else
    protocol.reconcenter = protocol.reconcenter(:)';
end
% windowcenter
if ~isfield(protocol, 'windowcenter')
    protocol.windowcenter = 0;
end
% windowwidth
if ~isfield(protocol, 'windowwidth')
    protocol.windowwidth = 100;
end
% reconkernel
if ~isfield(protocol, 'reconkernel')
    protocol.reconkernel = 'default';
end
% FOV
if ~isfield(protocol, 'reconFOV')
    % to fine FOV in pipe line
    tmp = findfield(pipe, 'FOV');
    if ~isempty(tmp)
        protocol.reconFOV = tmp;
    else
        % to find a default FOV?
        if isfield(system, 'maxFOV')
            protocol.reconFOV = system.maxFOV;
        else
            protocol.reconFOV = 0;
        end
    end
end
% Tube Angle
if ~isfield(protocol, 'TubeAngle')
    if isfield(protocol, 'startangle')
        % set TubeAngle with startangle
        protocol.TubeAngle = protocol.startangle;
    else
        protocol.TubeAngle = 0;
    end
end
% CollimatedSliceThickness
if ~isfield(protocol, 'CollimatedSliceThickness')
    protocol.CollimatedSliceThickness = CollimatedSliceThickness;
    % I know it was token from the protocol.collimator
end
% TotalCollimationWidth
if isfield(protocol, 'collimator')
    collitoken = regexp(protocol.collimator, '\<([\d \.]+)x', 'tokens');
else
    collitoken = [];
end
if isempty(collitoken)
    % what??
    SliceNumber = 0;
else
    SliceNumber = str2double(collitoken{1}{1});
end
if ~isfield(protocol, 'TotalCollimationWidth')
    protocol.TotalCollimationWidth = protocol.CollimatedSliceThickness*SliceNumber;
end
% TableFeedperRotation
if ~isfield(protocol, 'TableFeedperRotation')
    switch lower(protocol.scan)
        case 'axial'
            if isfield(protocol, 'shotcouchstep')
                protocol.TableFeedperRotation = abs(protocol.shotcouchstep);
            else
                protocol.TableFeedperRotation = 0;
            end
        case 'helical'
            if isfield(protocol, 'shotcouchstep') && isfield(protocol, 'rotationspeed')
                protocol.TableFeedperRotation = abs(protocol.couchspeed*protocol.rotationspeed);
            else
                protocol.TableFeedperRotation = 0;
            end
            
        otherwise
            protocol.TableFeedperRotation = 0;
    end
end
% ProtocolName
if ~isfield(protocol, 'ProtocolName')
    protocol.ProtocolName = 'QuickStart';
end
% ImageOrientationPatient
if ~isfield(protocol, 'ImageOrientationPatient')
    protocol.ImageOrientationPatient = [1 0 0  0 cos(protocol.gantrytilt*pi/180) sin(protocol.gantrytilt*pi/180)];
end
% PixelSpacing
if ~isfield(protocol, 'PixelSpacing')
    if ~isfield(protocol, 'reconFOV')
        error(['Error in missing recon paramter, at least one of these parameters shall be set: ''reconFOV'', ' ...
            '''FOV'' (in nodes) and ''PixelSpacing''']);
    end
    PixelSpacing = protocol.reconFOV/min(protocol.imagesize);
    protocol.PixelSpacing = [PixelSpacing PixelSpacing];
    % I know the image pixels must be square 
end

% datablock
if ~isfield(protocol, 'datablock')
    protocol.datablock = [];
end

% pipeline replicate
if ~isfield(protocol, 'pipelinereplicate')
    protocol.pipelinereplicate = [];
end

% pipeline nodes prepare
if ~isfield(protocol, 'pipelineprepare')
    % if to run the prepares before pipeline
    protocol.pipelineprepare = true;
end


end
