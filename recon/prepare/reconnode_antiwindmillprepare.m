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

% default configurable parameters
% down-sampling of the images in anti windmill
if ~isfield(prmflow.pipe.(nodename), 'downsample')
    % 2-times down-sampling
    prmflow.pipe.(nodename).downsample = 2;
end
% blur
if ~isfield(prmflow.pipe.(nodename), 'Gblur')
    prmflow.pipe.(nodename).Gblur = 0;
end
% TV 
if ~isfield(prmflow.pipe.(nodename), 'TVmu')
    prmflow.pipe.(nodename).TVmu = 0.03 + 0.03i;
end
if ~isfield(prmflow.pipe.(nodename), 'TVlambda')
    prmflow.pipe.(nodename).TVlambda = prmflow.pipe.(nodename).TVmu./3;
end
if ~isfield(prmflow.pipe.(nodename), 'TVCrange')
    prmflow.pipe.(nodename).TVCrange = [-inf inf] + [-inf inf].*1i;
end
if ~isfield(prmflow.pipe.(nodename), 'TVNiter')
    prmflow.pipe.(nodename).TVNiter = 40 + 40i;
end
if ~isfield(prmflow.pipe.(nodename), 'TVtol')
    prmflow.pipe.(nodename).TVtol = 0 + 0i;
end
if ~isfield(prmflow.pipe.(nodename), 'TVlogC')
    prmflow.pipe.(nodename).TVlogC = 4.0 + 4.0i;
end
if isfield(prmflow.pipe.(nodename), 'fixlimit')
    % no inf plz
    if any(isinf(prmflow.pipe.(nodename).fixlimit))
        prmflow.pipe.(nodename).fixlimit(isinf(prmflow.pipe.(nodename).fixlimit)) = 1e9;
    end
else
    prmflow.pipe.(nodename).fixlimit = 100 + 100i;
end
if ~isfield(prmflow.pipe.(nodename), 'fixsigma')
    prmflow.pipe.(nodename).fixsigma = 1e-3 + 1e-3i;
end
% if ~isfield(prmflow.pipe.(nodename), 'boundaryOpt')
%     prmflow.pipe.(nodename).boundaryOpt = 'none';
% end

% GPU
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    prmflow.pipe.(nodename) = everything2single(prmflow.pipe.(nodename), 'any', 'gpusingle');
else
    % to single
    prmflow.pipe.(nodename) = everything2single(prmflow.pipe.(nodename));
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe-line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    % the md is H-H.0.S or A.0.S 
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    prmflow.pipe.(nodename).pipeline.relystrategy = 'stingy';
    % viewrely is the boundary expandings for TV denoise
    if ~isfield(prmflow.pipe.(nodename), 'viewrely')
        prmflow.pipe.(nodename).viewrely = 16;
    end
    viewrely = prmflow.pipe.(nodename).viewrely;
    prmflow.pipe.(nodename).pipeline.viewrely = [viewrely viewrely];
    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
    % min input
    if ~isfield(prmflow.pipe.(nodename).pipeline, 'inputminlimit')
        prmflow.pipe.(nodename).pipeline.inputminlimit = 16;
    end

    % private buffer
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).imagebound = [];
end

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end 


