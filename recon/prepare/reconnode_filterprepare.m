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

% prm
nodeprm = prmflow.pipe.(status.nodename);
if isfield(nodeprm, 'upsampling')
    upsampling = nodeprm.upsampling;
else
    upsampling = false;
end
if isfield(nodeprm, 'upsampgamma') || upsampling
    if ~isfield(nodeprm, 'upsampgamma') || isempty(nodeprm.upsampgamma)
        prmflow.pipe.(status.nodename).upsampgamma = [0.7 0.8854];
    end
    upsampling = true;
end
prmflow.pipe.(status.nodename).upsampling = upsampling;

% upsampling
if upsampling
    prmflow.recon.upsampled = true;
    prmflow.recon.Npixel_up = prmflow.recon.Npixel*2;
    prmflow.recon.delta_d_up = prmflow.recon.delta_d/2;
    prmflow.recon.midchannel_up = round(prmflow.recon.midchannel*4-2)/2;

else
    prmflow.recon.upsampled = false;
    prmflow.recon.Npixel_up = prmflow.recon.Npixel;
    prmflow.recon.delta_d_up = prmflow.recon.delta_d;
    prmflow.recon.midchannel_up = prmflow.recon.midchannel;
end

% design filter
Npixel = prmflow.recon.Npixel_up;
midchannel = prmflow.recon.midchannel_up;
delta_d = prmflow.recon.delta_d_up;
prmflow.recon.filter.basicfilter = loadfilter(nodeprm, Npixel, delta_d);
prmflow.recon.filter.Hlen = length(prmflow.recon.filter.basicfilter);
prmflow.recon.filter.Npixel = Npixel;
prmflow.recon.filter.midchannel = midchannel;

% set filtered
prmflow.recon.filtered = true;

% pipe line
if pipeline_onoff
    % default
    prmflow.pipe.(nodename).pipeline = struct();
    % the filer is H.0.N or H.1.N
    if upsampling
        prmflow.pipe.(nodename).pipeline.kernellevel = 1;
        % ask datasize for next node
        prmflow.pipe.(nodename).pipeline.nextdatasize = double(prmflow.recon.Npixel_up * prmflow.recon.Nslice);
    else
        prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    end
    prmflow.pipe.(nodename).pipeline.relystrategy = 'none';
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end