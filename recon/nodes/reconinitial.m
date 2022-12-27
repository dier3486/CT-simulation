function [prmflow, status] = reconinitial(prmflow, status)
% recon initial
% [prmflow, status] = reconinitial(prmflow, status)

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
if ~iscell(status.reconcfg)
    reconcfg = status.reconcfg;
else
    reconcfg = status.reconcfg{status.seriesindex};
end

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

% clean
prmflow = iniprmclean(prmflow);

% ini GPU
status.GPUinfo = initialGPU(prmflow.system.GPUdeviceindex);

% series UID
status.seriesUID{status.seriesindex} = dicomuid;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function prmflow = iniprmclean(prmflow)
% to fill up the paramters which could be used in recon but not configured
% hard code

% collimatorexplain
if ~isfield(prmflow.system, 'collimatorexplain')
    prmflow.system.collimatorexplain = [];
elseif ischar(prmflow.system.collimatorexplain)
    prmflow.system.collimatorexplain = readcfgfile(prmflow.system.collimatorexplain);
else
    error('Illegal prmflow.system.collimatorexplain class: %s!', class(prmflow.system.collimatorexplain));
end
% filetagsrule
if ~isfield(prmflow.system, 'filetagsrule')
    prmflow.system.filetagsrule = struct();
elseif ischar(prmflow.system.filetagsrule)
    prmflow.system.filetagsrule = readcfgfile(prmflow.system.filetagsrule);
else
    error('Illegal prmflow.system.filetagsrule class: %s!', class(prmflow.system.filetagsrule));
end
% filematchrule
if ~isfield(prmflow.system, 'filematchrule')
    prmflow.system.filematchrule = struct();
elseif ischar(prmflow.system.filematchrule)
    prmflow.system.filematchrule = readcfgfile(prmflow.system.filematchrule);
else
    error('Illegal prmflow.system.filematchrule class: %s!', class(prmflow.system.filematchrule));
end
% GPU+
if ~isfield(prmflow.system, 'GPUdeviceindex')
    if gpuDeviceCount > 0
        prmflow.system.GPUdeviceindex = 1;
    else
        prmflow.system.GPUdeviceindex = 0;
    end
end
% set dicom dictionary
if isfield(prmflow.system, 'dicomdictionary') && ~isempty(prmflow.system.dicomdictionary)
    dicomdict('set', prmflow.system.dicomdictionary);
end

% protocol
prmflow.protocol.rawdata = prmflow.rawdata;     % copy rawdata path
prmflow.protocol = iniprotocolclean(prmflow.protocol, prmflow.pipe);

% IOstandard
if ~isfield(prmflow, 'IOstandard')
    prmflow.IOstandard = [];
end
% ini corrtable (always)
% if ~isfield(prmflow, 'corrtable')
%     prmflow.corrtable = struct();
% end
prmflow.corrtable = struct();
% ini recon (always)
prmflow.recon = struct();
% NOTE: sometimes we need to maintain data in prmflow for follow-up series, so we don't clean most of the informations in 
% prmflow when they already exist. But the prmflow.recon will always be cleaned.

end


function protocol = iniprotocolclean(protocol, pipe)
% to fill up the paramters in protocol
% hard code

% imagethickness
collitoken = regexp(protocol.collimator, 'x([\d \.]+)\>', 'tokens');
if isempty(collitoken)
    % what??
    CollimatedSliceThickness = 0;
else
    CollimatedSliceThickness = str2double(collitoken{1}{1});
end
if ~isfield(protocol, 'imagesize')
    % default imagesize 512x512
    protocol.imagesize = [512, 512];
elseif length(protocol.imagesize) == 1
    protocol.imagesize = [protocol.imagesize, protocol.imagesize];
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
    protocol.couchdirection = sign(sign(protocol.shotcouchstep) + sign(protocol.couchspeed));
end
% gantrytilt
if ~isfield(protocol, 'gantrytilt')
    protocol.gantrytilt = 0;
    % in 180 degree
end
% reconcenter
if ~isfield(protocol, 'reconcenter')
    protocol.reconcenter = [0 0];
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
    end
end
% Tube Angle
if ~isfield(protocol, 'TubeAngle')
    protocol.TubeAngle = protocol.startangle;
    % set TubeAngle with startangle
end
% CollimatedSliceThickness
if ~isfield(protocol, 'CollimatedSliceThickness')
    protocol.CollimatedSliceThickness = CollimatedSliceThickness;
    % I know it was token from the protocol.collimator
end
% TotalCollimationWidth
collitoken = regexp(protocol.collimator, '\<([\d \.]+)x', 'tokens');
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
            protocol.TableFeedperRotation = abs(protocol.shotcouchstep);
        case 'helical'
            protocol.TableFeedperRotation = abs(protocol.couchspeed*protocol.rotationspeed);
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
    PixelSpacing = protocol.reconFOV/min(protocol.imagesize);
    protocol.PixelSpacing = [PixelSpacing PixelSpacing];
    % I know the image pixels must be square 
end

end
