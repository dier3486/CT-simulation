function [dataflow, prmflow, status] = reconnode_antiringprepare(dataflow, prmflow, status)
% prepare node of antiring
% [dataflow, prmflow, status] = reconnode_antiringprepare(dataflow, prmflow, status);

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
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% recover the calibration table (which needs to do that was an unreleased version.)
if isfield(prmflow.corrtable, nodename)
    aringcorr = prmflow.corrtable.(nodename);
else
    aringcorr = struct();
end

% default values (we forbidened to set them in recon configure files)
if isfield(nodeprm, 'alpha')
    aringcorr.alpha = nodeprm.alpha;
elseif ~isfield(aringcorr, 'alpha')
    % anti-ring correction intensity
    aringcorr.alpha = 1.0;
    % default is 1
end
if isfield(nodeprm, 'Ntheta')
    aringcorr.Ntheta = nodeprm.Ntheta;
elseif ~isfield(aringcorr, 'Ntheta')
    % theta (samples on angle) number in r-theta transformation
    aringcorr.Ntheta = 192;
end
if isfield(nodeprm, 'rtheta_evenodd')
    aringcorr.rtheta_evenodd = nodeprm.rtheta_evenodd;
elseif ~isfield(aringcorr, 'rtheta_evenodd')
    % the flag of even/odd in r-theta transformation
    aringcorr.rtheta_evenodd = true;
    % default is true
end
if isfield(nodeprm, 'Crange')
    aringcorr.Crange = nodeprm.Crange;
elseif ~isfield(aringcorr, 'Crange')
    % the image value range to fix the ring artifacts
    aringcorr.Crange = single([-inf inf]);
end
if isfield(nodeprm, 'ringcut')
    aringcorr.ringcut = nodeprm.ringcut;
elseif ~isfield(aringcorr, 'ringcut')
    % the cut-off in the ring artifacts
    aringcorr.ringcut = single(20.0);
end
if isfield(nodeprm, 'ringmover')
    aringcorr.ringmover = nodeprm.ringmover;
elseif ~isfield(aringcorr, 'ringmover')
    % the filter width of meanmov and medianmov if ring finding
    aringcorr.ringmover = [5 5 5];  % int
end
if isfield(nodeprm, 'ringsection')
    aringcorr.ringsection = nodeprm.ringsection;
elseif ~isfield(aringcorr, 'ringsection')
    % the sections number in finding and fixing the ring artifacts
    aringcorr.ringsection = 2;  % int or single, [1, inf)
    % 1: anti whole ring, n: anti 1/n ring. The n need not to common with the Ntheta or Nviewprot, even non-integer n is
    % toleranced.
end
prmflow.corrtable.(nodename) = aringcorr;

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % the anti-ring correction is H.0.N
    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end 
end

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    % put corrtable to GPU
    prmflow.corrtable.(nodename) = putinGPU(prmflow.corrtable.(nodename));
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end