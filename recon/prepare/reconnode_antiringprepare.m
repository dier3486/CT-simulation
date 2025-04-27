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

% default settings
if ~isfield(nodeprm, 'alpha')
    % anti-ring correction intensity
    prmflow.pipe.(status.nodename).alpha = 1.0;
    % default is 1
end
if ~isfield(nodeprm, 'Ntheta')
    % theta (samples on angle) number in r-theta transformation
    prmflow.pipe.(status.nodename).Ntheta = 192;
end
if ~isfield(nodeprm, 'rtheta_evenodd')
    % the flag of even/odd in r-theta transformation
    prmflow.pipe.(status.nodename).rtheta_evenodd = true;
    % default is true
end
if ~isfield(nodeprm, 'Crange')
    % the image value range to fix the ring artifacts
    prmflow.pipe.(status.nodename).Crange = [-inf inf];
end
if ~isfield(nodeprm, 'Cnorm')
    % the normalization in 'restainting' the ring artifacts, some how should be the window-width in review the images.
    prmflow.pipe.(status.nodename).Cnorm = 200;
    % default is 200
end
if ~isfield(nodeprm, 'restcutoff')
    % the cut-off in 'restainting' the ring artifacts
    prmflow.pipe.(status.nodename).restcutoff = 0.1;
end
if ~isfield(nodeprm, 'ringfilter')
    % the filter in finding the ring artifacts
    prmflow.pipe.(status.nodename).ringfilter = [-0.5 1 -0.5];
end
if ~isfield(nodeprm, 'ringsection')
    % the sections number in finding and fixing the ring artifacts
    prmflow.pipe.(status.nodename).ringsection = 4;
end
if ~isfield(nodeprm, 'sectmethod')
    % the method to 'link' the sections
    prmflow.pipe.(status.nodename).sectmethod = 'spline';
    % linear or spline
end

% pipe line
if pipeline_onoff
    % pipeline console paramters, default
    prmflow.pipe.(nodename).pipeline = struct();
    % the anti-ring is H.0.N 
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end