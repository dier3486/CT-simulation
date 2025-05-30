function [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status)
% recon node, BP (after filter)
% [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% I know the prmflow.recon.method was filled by reconnode_BPprepare
switch lower(prmflow.recon.method)
    case 'axial2d'
        % 2D Axial (need update)
        [dataflow, prmflow, status] = reconnode_Axial2DBackprojection(dataflow, prmflow, status);
    case 'axial3d'
        % 3D Axial
        [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status);
    case {'helical', 'helical3d', 'helicalpiline'}
        % Helical and Helical-Piline
        [dataflow, prmflow, status] = reconnode_HelicalBackprojection(dataflow, prmflow, status);
    case {'cradle', 'cradlepiline'}
        % Cradle-Piline
        1;
        [dataflow, prmflow, status] = reconnode_CradleBackprojection(dataflow, prmflow, status);
    otherwise
        error('Sorry, the reconstruction %s is not supported yet!', prmflow.recon.method);
end


end