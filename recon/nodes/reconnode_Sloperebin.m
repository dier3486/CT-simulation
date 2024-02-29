function [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status)
% recon node, Axial 'slope' rebin 
% [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);
% Support (X)DFS, gantry tilt, NO QDO

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

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_sloperebinprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
if isfield(nodeprm, 'pipeline_onoff')
    pipeline_onoff = status.pipeline_onoff & nodeprm.pipeline_onoff;
else
    pipeline_onoff = status.pipeline_onoff;
end

%1 slope fan-Radial

% pipeline consol
% read the data from input pool
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    % datasize
    datasize_toread = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % check the buffer left
    outpoolsize = dataflow.buffer.(nodename).outpoolsize - dataflow.buffer.(nodename).WritePoint + 1;
    datasize_toread = min(datasize_toread, outpoolsize);

    % copy pool to buffer.outpool
    [dataflow.buffer.(nodename).outpool, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow.buffer.(nodename).outpool, statuscurr.ReadPoint, dataflow.buffer.(nodename).WritePoint, ...
        datasize_toread, {}, true);

    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + datasize_toread; 
    % move the buffer.outpool's write point
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint + datasize_toread;
    % Note: the buffer.(nodename).outpool is a private buffer, not the next pool. The data in buffer.(nodename).outpool will be
    % copied to the next pool after the correcion works done.
end

[dataflow, prmflow, ~] = reconnode_SlopeRadialrebin(dataflow, prmflow, status);

%2 slope Azi
[dataflow, prmflow, ~] = reconnode_SlopeAzirebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end