function [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status)
% recon node, off-focal correction
% [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status);
% only for axial

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

scantype = prmflow.raw.scan;
switch lower(scantype)
    case {'axial', 'static'}
%         [dataflow, prmflow, status] = reconnode_offfocalaxial(dataflow, prmflow, status);
        [dataflow, prmflow, status] = reconnode_offfocalcorr_old(dataflow, prmflow, status);
        
    case {'helical', 'halfaxial'}
        [dataflow, prmflow, status] = reconnode_offfocalhelical(dataflow, prmflow, status);
    otherwise
        % what?
        warning('Illeagal scan type for off-focal correction: %s!', scantype);
end

% status (tmp)
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end