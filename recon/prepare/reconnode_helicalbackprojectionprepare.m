function [dataflow, prmflow, status] = reconnode_helicalbackprojectionprepare(dataflow, prmflow, status)
% prepare node, just the reconnode_backprojectionprepare
% [dataflow, prmflow, status] = reconnode_helicalbackprojectionprepare(dataflow, prmflow, status);

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

if ~isfield(prmflow.recon, 'method') || isempty(prmflow.recon.method)
    prmflow.recon.method = 'helical';
end
[dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status);

end