function [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status)
% recon node, nonlinear correction
% [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status);
% the algorithm is exactly same as beamharden corr, so we suggest to assign beamharden and nonlinear correction in one node in
% nodesentry.m

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
Nview = prmflow.raw.Nview;
% Nfocal = prmflow.raw.Nfocal;

% calibration table
nonlcorr = prmflow.corrtable.(status.nodename);
nonlorder = nonlcorr.order(1);
nonlpoly = reshape(nonlcorr.main, nonlcorr.Npixel*nonlcorr.Nslice, nonlorder, []);
% DFS
if isfield(nonlcorr, 'focalnumber') && nonlcorr.focalnumber
    Nfocal = nonlcorr.focalnumber;
else
    Nfocal = size(nonlpoly, 3);
end

% beam harden polynomial
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
for ifocal = 1:Nfocal
    dataflow.rawdata(:, ifocal:Nfocal:end) = iterpolyval(nonlpoly(:, :, ifocal), dataflow.rawdata(:, ifocal:Nfocal:end));
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end