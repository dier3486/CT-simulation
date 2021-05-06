function [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status)
% recon node, log2
% [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);

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

% Z0
if isfield(prmflow, 'system') && isfield(prmflow.system, 'DBBzero')
    Z0 = prmflow.system.DBBzero;
else
    Z0 = 16384;
end
% offset
if isfield(dataflow, 'offset')
    dataflow.rawdata = dataflow.rawdata - mean(dataflow.offset.rawdata, 2);
else
    dataflow.rawdata = dataflow.rawdata - Z0;
end
dataflow.rawdata(dataflow.rawdata<=0) = nan;
% log2
dataflow.rawdata = -log2(dataflow.rawdata) + log2(single(dataflow.rawhead.Integration_Time));

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end