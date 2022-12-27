function [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status)
% recon node, Axial rebin 
% [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
% just put the reconnode_rebinprepare, reconnode_Azirebin and
% reconnode_Radialrebin in one function.
% Support QDO, (X)DFS and QDO+DFS
% mainly use in 2D recon, for 3D recon try the Sloperebin could be better

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

% % test DFS
% prmflow.recon.focalspot = [2 3];
% prmflow.recon.Nfocal = 2;
1;
% rebin prepare
[prmflow, ~] = reconnode_axialrebinprepare(prmflow, status);

% Azi rebin
[dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);

% Radial rebin
[dataflow, prmflow, ~] = reconnode_Radialrebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end