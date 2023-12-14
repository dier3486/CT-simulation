function [prmflow, status] = reconnode_upsampleprepare(prmflow, status)
% prepare node, X-upsample prepare
% [prmflow, status] = reconnode_upsampleprepare(prmflow, status);

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

prmflow.recon.upsampled = true;

prmflow.recon.Npixel_up = prmflow.rebin.Nreb*2;
prmflow.recon.delta_d_up = prmflow.rebin.delta_d/2;
prmflow.recon.midchannel_up = round(prmflow.rebin.midU_phi*4-2)/2;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end