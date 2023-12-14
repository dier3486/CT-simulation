function [dataflow, prmflow, status] = reconnode_filterprepare(dataflow, prmflow, status)
% prepare node, filter prepare
% [prmflow, status] = reconnode_filterprepare(prmflow, status);

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

% nodename
nodename = status.nodename;
% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if isfield(prmflow.recon, 'upsampled') && prmflow.recon.upsampled
    % has been upsampled
    Npixel = prmflow.recon.Npixel_up;
    midchannel = prmflow.recon.midchannel_up;
    delta_d = prmflow.recon.delta_d_up;
else
    % not updampled
    Npixel = prmflow.rebin.Nreb;
    midchannel = prmflow.rebin.midU_phi;
    delta_d = prmflow.rebin.delta_d;
end

% the filter is
filter = prmflow.pipe.(status.nodename);

% design filter
prmflow.recon.filter.basicfilter = loadfilter(filter, Npixel, delta_d);
prmflow.recon.filter.Hlen = length(prmflow.recon.filter.basicfilter);
prmflow.recon.filter.Npixel = Npixel;
prmflow.recon.filter.midchannel = midchannel;

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    dataflow.buffer.(nodename) = struct();
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end