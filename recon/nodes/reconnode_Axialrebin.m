function [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status)
% recon node, Axial-2D rebin 
% [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
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

% rebin prepare
% [prmflow, ~] = reconnode_axialrebinprepare(prmflow, status);
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_axialrebinprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

[dataflow, prmflow, status] = reconnode_Rebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end