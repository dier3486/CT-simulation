function [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status)
% recon node, interp the rawdata to equal-fanangle, 
% [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status);

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
    [prmflow, status] = reconnode_fanangleprepare(prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters to use in prmflow
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
Nfocal = prmflow.raw.Nfocal;

equalfaninterp = prmflow.rebin.equalfaninterp;
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nfocal, []);

% interp1
for ii = 1:Nslice*Nfocal
    dataflow.rawdata(:, ii, :) = interp1(dataflow.rawdata(:, ii, :), equalfaninterp(:, ii));
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end