function [dataflow, prmflow, status] = reconnode_datareloadprepare(dataflow, prmflow, status)
% prepare node, datareload prepare
% [dataflow, prmflow, status] = reconnode_datareloadprepare(dataflow, prmflow, status);

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

if ~status.pipeline.(status.nodename).pipeline_onoff
    status.pipeline.(status.nodename).pipeline_onoff = true;
    warning('The node datareload can only be used in pipeline-on, whose pipeline_onoff is forcedly set to true!');
end

% nodename
nodename = status.nodename;

% datablock (size)
if ~isfield(prmflow.pipe.(nodename), 'datablock') || isempty(prmflow.pipe.(nodename).datablock)
    prmflow.pipe.(nodename).datablock = 512;
end
% data fields
if isfield(prmflow.pipe.(nodename), 'reloadfields')
    if ~iscell(prmflow.pipe.(nodename).reloadfields)
        prmflow.pipe.(nodename).reloadfields = regexp(prmflow.pipe.(nodename).reloadfields, '(, +)|(,)', 'split');
    end
else
    % default reloadfields are rawdata and rawhead
    prmflow.pipe.(nodename).reloadfields = {'rawdata', 'rawhead'};
end
% But only the rawdata and rawhead can be reload now. (image is not supported yet)
% clear afer reloading?
if ~isfield(prmflow.pipe.(nodename), 'clearafterreload') || isempty(prmflow.pipe.(nodename).clearafterreload)
    prmflow.pipe.(nodename).clearafterreload = true;
    % the reloadfields will be deleted after reloading
end

% input pool
dataflow.pipepool.(nodename) = rmfield(status.defaultpool, 'data');
dataflow.pipepool.(nodename).datafields = prmflow.pipe.(nodename).reloadfields;
dataflow.pipepool.(nodename).poolsize = 0;  % writing forbidden
dataflow.pipepool.(nodename).recylestrategy = 0;
% The input pool is used to record the points in reloading data. What the data will be reloading is saved in dataflow but not
% the inputpool.

% private buffer
dataflow.buffer.(nodename) = struct();
dataflow.buffer.(nodename).iblock = 1;
dataflow.buffer.(nodename).ishot = 0;
if isfield(prmflow, 'prepare')
    dataflow.buffer.(nodename).Nview = prmflow.prepare.Nview;
    dataflow.buffer.(nodename).viewpershot = prmflow.prepare.viewpershot;
    dataflow.buffer.(nodename).Nshot = prmflow.prepare.Nshot;
    dataflow.buffer.(nodename).startshot = prmflow.prepare.startshot;
else
    dataflow.buffer.(nodename).Nview = inf;
    dataflow.buffer.(nodename).viewpershot = inf;
    dataflow.buffer.(nodename).Nshot = 1;
    dataflow.buffer.(nodename).startshot = 1;
end
% Note: this node is depending on the prmflow.prepare to get the shot related parameters, therefore the previous nodes shall
% keep the prmflow.prepare in updating to support the datareload, or it could lay in error.

% node.pipeline
prmflow.pipe.(nodename).pipeline.kernellevel = 1;    % do not copy curr to next
prmflow.pipe.(nodename).pipeline.inputminlimit = 0;  % need not input data
prmflow.pipe.(nodename).pipeline.priostep = false;

end