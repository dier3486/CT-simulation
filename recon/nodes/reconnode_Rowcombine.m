function [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status)
% recon node, row combine by detector.rowcombine, in working
% [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status);

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
% Nview = prmflow.raw.Nview;
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;

if isfield(prmflow.system.detector, 'rowcombine') && ~isempty(prmflow.system.detector, 'rowcombine')
    rowcombine = prmflow.system.detector.rowcombine;
    % combine the rawdata
    [dataflow.rawdata, Nrowcombined] = detectorslicemerge(dataflow.rawdata, Npixel, Nslice, rowcombine, 'mean');
    % combine the detector position
    [prmflow.system.detector.position, ~] = ...
        detectorslicemerge(prmflow.system.detector.position, Npixel, Nslice, rowcombine, 'mean');
    prmflow.raw.Nslice = Nrowcombined;
    % shall in prepare
    prmflow.rebin.Nslice = Nrowcombined;
    prmflow.recon.Nslice = Nrowcombined;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end