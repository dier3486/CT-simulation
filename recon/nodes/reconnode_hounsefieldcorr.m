function [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status)
% recon node, Hounsefield Units correction
% [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status);

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

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
% parameters set in pipe
HCprm = prmflow.pipe.(status.nodename);

if isfield(HCprm, 'HCscale')
    HCscale = HCprm.HCscale;
else
    HCscale = 1000;
end

% scale
dataflow.rawdata = dataflow.rawdata.*HCscale;

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    HUcorr = prmflow.corrtable.(status.nodename);
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview).*HUcorr.main(:)';
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end