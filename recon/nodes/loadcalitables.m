function [prmflow, status] = loadcalitables(prmflow, status)
% load calibration tables
% [prmflow, status] = loadcalitables(reconcfg, prmflow, status);

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


% corrpath is the path to looking for corr files
if isfield(prmflow, 'corrpath')  && ~isempty(prmflow.corrpath)
    corrpath = prmflow.corrpath;
elseif isfield(prmflow.system, 'corrpath') && ~isempty(prmflow.system.corrpath)
    corrpath = prmflow.system.corrpath;
else
    [corrpath, ~, ~] = fileparts(prmflow.rawdata);
end
% corrext is the ext of the corr files
if isfield(prmflow.system, 'corrext')
    corrext = prmflow.system.corrext;
    % NOTE: set corrext = '.(corr|ct)' to looking for both .corr and .ct files
else
    % default corrext is .corr
    corrext = '.corr';
end

% detector position
% try to look for detector corr in corrpath
detcorrfile = corrcouplerule(prmflow.protocol, corrpath, prmflow.system.filematchrule, 'detector', corrext);
if isempty(detcorrfile)
    if isfield(prmflow, 'detector_corr')
        detcorrfile = prmflow.detector_corr;
    else
        detcorrfile = prmflow.system.detector_corr;
    end
end
prmflow.system.detector_corr = detcorrfile;
% load detector
prmflow.system.detector = loaddetector(detcorrfile, prmflow.IOstandard, prmflow.protocol.collimator, ...
    prmflow.system.collimatorexplain.collimator);

% put focalposition in system
if ~isfield(prmflow.system, 'focalposition') || isempty(prmflow.system.focalposition)
    prmflow.system.focalposition = prmflow.system.detector.focalposition;
end

% to prm.recon
prmflow.recon.Nslice = prmflow.system.detector.Nmergedslice;
prmflow.recon.Npixel = double(prmflow.system.detector.Npixel);

% other tables
pipenodes = fieldnames(prmflow.pipe);
for ii = 1:length(pipenodes)
    % node name
    nodename = regexp(pipenodes{ii}, '_', 'split');
    nodename = lower(nodename{1});
    % is the .corr has been defined?
    if ~isfield(prmflow.pipe.(pipenodes{ii}), 'corr') || isempty(prmflow.pipe.(pipenodes{ii}).corr)
        % if the nodename is included in filematchrule to find the corr file
        if isfield(prmflow.system.filematchrule, nodename)
            corrfile = corrcouplerule(prmflow.protocol, corrpath, prmflow.system.filematchrule, nodename, corrext);
            if corrfile
                prmflow.pipe.(pipenodes{ii}).corr = corrfile;
%                 corrfile
            end
        end
    end
    % load the corrfile
    if isfield(prmflow.pipe.(pipenodes{ii}), 'corr')
        if isempty(prmflow.pipe.(pipenodes{ii}).corr)
            warning('Not found the calibration table of recon node: %s!', pipenodes{ii});
            continue;
        end
        prmflow.corrtable.(pipenodes{ii}) = loaddata(prmflow.pipe.(pipenodes{ii}).corr, prmflow.IOstandard, nodename);
        prmflow.corrtable.(pipenodes{ii}).filename = prmflow.pipe.(pipenodes{ii}).corr;
        % reuse corr for different collimator
        prmflow.corrtable.(pipenodes{ii}) = ...
            collimatedcorr(prmflow.corrtable.(pipenodes{ii}), nodename, prmflow.system.detector);
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function detector = loaddetector(detcorrfile, IOstandard, collimator, collimatorexplain)

% load corr file (detector)
det_corr = loaddata(detcorrfile, IOstandard, 'detector');
% I know the arg 'detector' is for CRIS .ct files. 
% explain the collimator
detector = collimatorexposure(collimator, [], det_corr, collimatorexplain);
% mergeslice
detector.position = reshape(detector.position, [], 3);
[detector.position, detector.Nmergedslice] = ...
    detectorslicemerge(detector.position, detector.Npixel, detector.Nslice, detector.slicemerge, 'mean');
% copy other parameters from det_corr
detector = structmerge(detector, det_corr);
% reshape focal position
if isfield(detector, 'focalposition')
    detector.focalposition = reshape(detector.focalposition, [], 3);
end
% pixel range
if isfield(detector, 'pixelrange')
    detector.pixelrange = reshape(detector.pixelrange, 2, []);
    detector.Nprange = max(mod(detector.pixelrange(2, :)-detector.pixelrange(1, :), detector.Npixel) + 1);
end

end

