function [dataflow, prmflow, status] = reconnode_Helicalrebin(dataflow, prmflow, status)
% recon node, Helical rebin 
% [dataflow, prmflow, status] = reconnode_Helicalrebin(dataflow, prmflow, status);
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
    [dataflow, prmflow, status] = reconnode_helicalrebinprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% call Rebin node
[dataflow, prmflow, status] = reconnode_Rebin(dataflow, prmflow, status);

end
