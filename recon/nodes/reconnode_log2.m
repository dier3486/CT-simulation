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

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% offset
if ~isfield(prmflow.correction, 'offset')
    % save the offset to prmflow.raw.offset
    if isfield(dataflow, 'offset') && isfield(dataflow.offset, 'rawdata')
        prmflow.correction.offset = mean(dataflow.offset.rawdata, 2);
        % Note: the offset will not be put in dataflow.pipepool.(nodename)!
    else
        if isfield(prmflow, 'system') && isfield(prmflow.system, 'DBBzero')
            prmflow.correction.offset = single(prmflow.system.DBBzero);
        else
            prmflow.correction.offset = single(16384);
        end
    end
end

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
end

% main
log2Kernelfuntion();

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function log2Kernelfuntion()
        % The anonymous function is static
        debug = [];

        if pipeline_onoff
            nextnode = status.currentjob.nextnode;
            carrynode = status.currentjob.carrynode;
            if isempty(nextnode) || strcmpi(nextnode, 'NULL')
                return;
            end
        end

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(dataflow.pipepool.(nextnode), plconsol.Index_out);
%             rawdata = dataflow.pipepool.(carrynode).data.rawdata(:, index_out);
            Integration_Time = single(dataflow.pipepool.(carrynode).data.rawhead.Integration_Time(index_out));
        else
            index_out = 1:prmflow.raw.Nview;
%             rawdata = dataflow.rawdata;
            Integration_Time = single(dataflow.rawhead.Integration_Time);
        end

        if pipeline_onoff
            dataflow.pipepool.(carrynode).data.rawdata(:, index_out) = ...
                log2opterators(dataflow.pipepool.(carrynode).data.rawdata(:, index_out), prmflow.correction.offset, Integration_Time);
        else
            dataflow.rawdata = log2opterators(dataflow.rawdata, prmflow.correction.offset, Integration_Time);
            % jobdone
            status.jobdone = true;
        end

    end
end


function rawdata = log2opterators(rawdata, offset, Integration_Time)
% offset
rawdata = rawdata - offset;

% negative value
Sneg = rawdata<=0;
% if any(Sneg(:))
rawneg = -rawdata(Sneg)./log(2) + 1/log(2);
1;
% log2
rawdata = -log2(rawdata);
rawdata(Sneg) = rawneg;

% intigration time
rawdata = rawdata + log2(Integration_Time);
end
