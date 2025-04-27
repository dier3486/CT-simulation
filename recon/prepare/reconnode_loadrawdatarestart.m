function [dataflow, prmflow, status] = reconnode_loadrawdatarestart(dataflow, prmflow, status)
% read rawdata restart (for 'infinite' pipe-line series)
% [dataflow, prmflow, status] = reconnode_loadrawdatarestart(dataflow, prmflow, status)

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

% left view number
prmflow.raw.viewleft = prmflow.raw.Nview - prmflow.raw.viewread;

% update viewnumber and Nview
prmflow.raw.viewnumber = prmflow.raw.viewnumber - prmflow.raw.viewread;
prmflow.raw.Nview = min(prmflow.raw.viewnumber, prmflow.raw.maxviewnumber);
% reset viewread
prmflow.raw.viewread = 0;

% blocked reading
if prmflow.raw.datablock_onoff
    datablock = prmflow.protocol.datablock;
    Nblock = ceil(prmflow.raw.Nview / datablock);
    prmflow.raw.Nblock = Nblock;
    prmflow.raw.datablocksize = ...
        [repmat(datablock, 1, Nblock-1)  prmflow.raw.Nview - datablock * (Nblock-1)];
else
    % no blocks
    prmflow.raw.Nblock = 1;
    prmflow.raw.datablocksize = prmflow.raw.Nview;
end

if isinf(prmflow.protocol.shotnumber)
    % TBC
end

% ini iblock
prmflow.raw.iblock = 1;

% job done
status.jobdone = true;

end