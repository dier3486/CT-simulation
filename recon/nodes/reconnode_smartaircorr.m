function [dataflow, prmflow, status] = reconnode_smartaircorr(dataflow, prmflow, status)
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
    [dataflow, prmflow, status] = reconnode_smartairprepare(dataflow, prmflow, status);
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
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
end

% main
if pipeline_onoff
    dataflow.pipepool.(carrynode).data = ...
        smartaircorrKernelfuntion(dataflow.pipepool.(nextnode), dataflow.pipepool.(carrynode).data);
else
    dataflow = smartaircorrKernelfuntion([], dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function data = smartaircorrKernelfuntion(nextpool, data)
        % The anonymous function is static
        debug = [];
        
        % prm
        anglepercode = prmflow.raw.anglepercode;
        angulationzero = prmflow.raw.angulationzero;
        % calibration table
        smartaircorr = prmflow.corrtable.(status.nodename);
        % GPU on/off
        GPUonoff = status.currentjob.GPUdevice > 0;

        % data classes
        if GPUonoff
            dataclass_single = ones(1, 'single', 'gpuArray');
%             dataclass_raw = gpuArray(ones(1, 'like', dataIn.rawdata));
        else
            dataclass_single = ones(1, 'single');
%             dataclass_raw = ones(1, 'like', dataIn.rawdata);
        end

        % pipeline consol
        if pipeline_onoff
            % index
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(nextpool, plconsol.Index_out);
            Nview = length(index_out);
        else
            Nview = size(data.rawdata, 2);
            index_out = 1:Nview;
        end

        % KVmA
        % KVmA = cast(data.rawhead.KV(index_out).* data.rawhead.mA(index_out), 'like', dataclass_single);
        % viewangle
        viewangle = cast(double(data.rawhead.AngleEncoder(index_out)).*anglepercode + angulationzero, ...
            'like', dataclass_single);

        % air corr
        data.rawdata(:, index_out) = smartaircorrwithoutref(data.rawdata(:, index_out), prmflow.raw, viewangle, smartaircorr);
        
        if ~pipeline_onoff
            % jobdone
            status.jobdone = true;
        end

    end

end
