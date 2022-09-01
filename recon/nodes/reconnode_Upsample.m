function [dataflow, prmflow, status] = reconnode_Upsample(dataflow, prmflow, status)
% recon node, double upsample
% [dataflow, prmflow, status] = reconnode_Upsample(dataflow, prmflow, status);

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

% prm
nodeprm = prmflow.pipe.(status.nodename);
if isfield(nodeprm, 'upsampgamma') && ~isempty(nodeprm.upsampgamma)
    upsampgamma = nodeprm.upsampgamma;
else
    upsampgamma = [0.7 0.8854];
end
Npixel = prmflow.recon.Npixel;

% double up
dataflow.rawdata = doubleup(reshape(dataflow.rawdata, Npixel, []), upsampgamma);

% up prm
prmflow.recon.Npixel = Npixel*2;
prmflow.recon.midchannel = prmflow.recon.midchannel*2-1;
prmflow.recon.delta_d = prmflow.recon.delta_d/2;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end