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

%1 slope fan-Radial
[dataflow, prmflow, ~] = reconnode_SlopeRadialrebin(dataflow, prmflow, status);

%2 slope Azi
[dataflow, prmflow, ~] = reconnode_SlopeAzirebin(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end