function [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status);
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
aircorrKernelfuntion();

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function aircorrKernelfuntion()
        % The anonymous function is static
        debug = [];
        
        if pipeline_onoff
            nextnode = status.currentjob.nextnode;
            carrynode = status.currentjob.carrynode;
            if isempty(nextnode) || strcmpi(nextnode, 'NULL')
                return;
            end
        end

        % calibration table
        aircorr = prmflow.corrtable.(status.nodename);
        % refnumber
        refnumber = prmflow.correction.air.refnumber;

        if pipeline_onoff
            % index
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(dataflow.pipepool.(nextnode), plconsol.Index_out);
            % KVmA
            KVmA = dataflow.pipepool.(carrynode).data.rawhead.KV(index_out).* ...
                dataflow.pipepool.(carrynode).data.rawhead.mA(index_out);
            viewangle = dataflow.pipepool.(carrynode).data.rawhead.viewangle(index_out);

            % air corr
            dataflow.pipepool.(carrynode).data.rawdata(:, index_out) = ...
                aircorrwithoutref(dataflow.pipepool.(carrynode).data.rawdata(:, index_out), prmflow, KVmA, viewangle, aircorr);

            % ini refblock
            if ~isfield(dataflow.pipepool.(carrynode).data.rawhead, 'refblock')
                dataflow.pipepool.(carrynode).data.rawhead.refblock = ...
                    false(refnumber, size(dataflow.pipepool.(carrynode).data.rawdata, 2));
            end
        else
            % KVmA
            KVmA = dataflow.rawhead.KV.*dataflow.rawhead.mA;
            viewangle = dataflow.rawhead.viewangle;

            % air corr
            dataflow.rawdata = aircorrwithoutref(dataflow.rawdata, prmflow, KVmA, viewangle, aircorr);

            % ini refblock
            dataflow.rawhead.refblock = false(refnumber, size(dataflow.rawdata, 2));

            % jobdone
            status.jobdone = true;
        end

    end

end



