function [dataflow, prmflow, status] = reconnode_backprojectionrestart(dataflow, prmflow, status)
% BP restart
% [dataflow, prmflow, status] = reconnode_backprojectionrestart(dataflow, prmflow, status);

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

% switch recon method
switch lower(prmflow.recon.method)
    case {'axial2d', 'axial3d'}
        % TBC
    case {'helical', 'helical3d', 'helicalpiline'}
        % helical is always 3D
        [dataflow, prmflow, status] = reconnode_helicalbackprojectionrestart(dataflow, prmflow, status);
    case {'axialhalf'}
        0;
        % TBC
    otherwise
        % do nothing
        1;
        % no topo
end

end
