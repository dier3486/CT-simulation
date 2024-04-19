function [dataflow, prmflow, status] = reconnode_materialdecompprepare(dataflow, prmflow, status)
% prepare node, two-material decomposition prepare
% [dataflow, prmflow, status] = reconnode_materialdecompprepare(dataflow, prmflow, status);

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

if isfield(prmflow.system, 'slicezebra')
    slicezebra = prmflow.system.slicezebra;
else
    slicezebra = false;
end
if ~slicezebra
    error('Only slice zebra mode now!');
end

% recover the calibration table (which needs to do that was an unreleased version.)
% calibration table to single
mdcorr = everything2single(prmflow.corrtable.(status.nodename));

% tableZr
if ~isfield(mdcorr, 'tableZr')
    mdcorr.tableZr = (mdcorr.Tablelk2r.*(mdcorr.ZB^3 - mdcorr.ZA^3)+mdcorr.ZA^3).^(1/3);
end

% copy back
prmflow.corrtable.(status.nodename) = mdcorr;

% create prmflow.materialdecomp for posterior material-maps
prmflow.materialdecomp.slicezebra = slicezebra;
prmflow.materialdecomp.HU = 1000;
prmflow.materialdecomp.ZA = mdcorr.ZA;
prmflow.materialdecomp.ZB = mdcorr.ZB;
prmflow.materialdecomp.densityA = mdcorr.densityA;
prmflow.materialdecomp.densityB = mdcorr.densityB;
prmflow.materialdecomp.elecdensity_flag = mdcorr.elecdensity_flag;

status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end 