function [dataflow, prmflow, status] = reconnode_antiwindmillprepare(dataflow, prmflow, status)
% prepare node, anti windmill artifact prepare
% [dataflow, prmflow, status] = reconnode_antiwindmillprepare(dataflow, prmflow, status);

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

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% recover the calibration table (which needs to do that was an unreleased version.)
if isfield(prmflow.corrtable, nodename)
    awindcorr = prmflow.corrtable.(nodename);
else
    awindcorr = struct();
end

% default values (we forbidened to set them in recon configure files)
% down-sampling of the images in anti windmill
if ~isfield(awindcorr, 'downsample')
    % 2-times down-sampling
    awindcorr.downsample = 2;
end
% blur
if ~isfield(awindcorr, 'Gblur')
    awindcorr.Gblur = single(0);
end
% TV 
if ~isfield(awindcorr, 'TVmu')
    awindcorr.TVmu = single(0.03 + 0.03i);
end
if ~isfield(awindcorr, 'TVlambda')
    awindcorr.TVlambda = awindcorr.TVmu./3;
end
if ~isfield(awindcorr, 'TVCrange')
    awindcorr.TVCrange = single([-inf inf] + [-inf inf].*1i);
end
if ~isfield(awindcorr, 'TVNiter')
    awindcorr.TVNiter = single(40 + 40i);
end
if ~isfield(awindcorr, 'TVtol')
    awindcorr.TVtol = single(0 + 0i);
end
if ~isfield(awindcorr, 'TVlogC')
    awindcorr.TVlogC = single(4.0 + 4.0i);
end
if isfield(awindcorr, 'fixlimit')
    % no inf plz
    if any(isinf(awindcorr.fixlimit))
        awindcorr.fixlimit(isinf(awindcorr.fixlimit)) = single(1e9);
    end
else
    awindcorr.fixlimit = single(100 + 100i);
end
if ~isfield(awindcorr, 'fixsigma')
    awindcorr.fixsigma = single(1e-3 + 1e-3i);
end
if ~isfield(awindcorr, 'relybuffer')
    awindcorr.imagerely = 16;
end
% if ~isfield(awindcorr, 'boundaryOpt')
%     awindcorr.boundaryOpt = 'none';
% end
prmflow.corrtable.(nodename) = awindcorr;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe-line
if pipeline_onoff
    % the aw is H-H.0.S
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    prmflow.pipe.(nodename).pipeline.relystrategy = 'stingy';
    % viewrely(imagerely) is the boundary expandings for TV denoise
    if isfield(prmflow.pipe.(nodename), 'imagerely')
        viewrely = prmflow.pipe.(nodename).imagerely;
    else
        viewrely = prmflow.corrtable.(nodename).imagerely;
    end
    prmflow.pipe.(nodename).pipeline.viewrely = [viewrely viewrely];
    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
    % min input
    if ~isfield(prmflow.pipe.(nodename).pipeline, 'inputminlimit')
        prmflow.pipe.(nodename).pipeline.inputminlimit = 16;
    end

    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end 

    % private buffer
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).imagebound = [];
end

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    % put corrtable to GPU
    prmflow.corrtable.(nodename) = putinGPU(prmflow.corrtable.(nodename));
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end 


