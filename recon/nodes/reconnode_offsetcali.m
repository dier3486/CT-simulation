function [dataflow, prmflow, status] = reconnode_offsetcali(dataflow, prmflow, status)
% recon node, offset cali
% [dataflow, prmflow, status] = reconnode_offsetcali(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_offsetcaliprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

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
offsetcaliKernelfuntion();

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function offsetcaliKernelfuntion()
        % The anonymous function is static
        debug = [];
        
        nodename = status.nodename;
        nodeprm = prmflow.pipe.(nodename);
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
            isshotend = plconsol.isshotend;
        else
            isshotend = true;
        end
        
        % main in buffer
        if pipeline_onoff
            dataflow.buffer.(nodename) = ...
                offsetcali(dataflow.buffer.(nodename), dataflow.pipepool.(carrynode).data.rawdata(:, index_out));
        else
            dataflow.buffer.(nodename) = ...
                offsetcali(dataflow.buffer.(nodename), dataflow.rawdata);
        end
        
        % create offcorr
        if isshotend
            % paramters for offcorr
            dataflow.calibration.log2corr = caliprmforcorr(prmflow, nodeprm.corrversion);
            dataflow.calibration.log2corr.main = dataflow.buffer.(nodename).main;
            dataflow.calibration.log2corr.variance = dataflow.buffer.(nodename).variance;
            dataflow.calibration.log2corr.Nslice = prmflow.raw.Nslice;
            dataflow.calibration.log2corr.mainsize = length(dataflow.buffer.(nodename).main(:));
            dataflow.calibration.log2corr.reserved = [];
        end

        if ~pipeline_onoff
            % jobdone
            status.jobdone = true;
        end

    end
end
