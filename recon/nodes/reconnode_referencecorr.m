function [dataflow, prmflow, status] = reconnode_referencecorr(dataflow, prmflow, status)
% recon node, (air) reference correction
% [dataflow, prmflow, status] = reconnode_referencecorr(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% Axial or helical
switch lower(prmflow.raw.scan)
    case 'axial'
        [dataflow, prmflow, status] = reconnode_referenceaxialcorr(dataflow, prmflow, status);
    case {'helical', 'static', 'halfaxial'}
        [dataflow, prmflow, status] = reconnode_referencehelicalcorr(dataflow, prmflow, status);
    otherwise
        status.jobdone = false;
        status.errorcode = 1;
        status.errormsg = sprintf('Not support scan mode %s in air correcion', prmflow.raw.scan);
end

end
