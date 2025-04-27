function [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status)
% recon node, nonlinear correction (beamharden correction)
% [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status);
% The beamharden correction and nonlinear correction is a same node in nodesentry.m. They have different name different
% calibration table but exactly same algorithim in recon correction.

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

% no prepare of non-linear (and beamharden) correction

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
if pipeline_onoff
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
    dataflow.pipepool.(carrynode).data = ...
        nonlinearcorrKernelfuntion(dataflow.pipepool.(nextnode), dataflow.pipepool.(carrynode).data);
else
    dataflow = nonlinearcorrKernelfuntion([], dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function data = nonlinearcorrKernelfuntion(nextpool, data)
        % The anonymous function is static
        debug = [];

        % calibration table
        nonlcorr = prmflow.corrtable.(status.nodename);
        nonlorder = nonlcorr.order(1);
        
        % DFS
        Nfocal = prmflow.raw.Nfocal;        
        nonlpoly = reshape(nonlcorr.main, [], nonlorder, nonlcorr.focalnumber);
        
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

        % beam harden polynomial
        for ifocal = 1:Nfocal
            viewindex = index_out(ifocal : Nfocal : Nview);
            ifocal_mod = mod((ifocal-1), nonlcorr.focalnumber) + 1;
            data.rawdata(:, viewindex) = iterpolyval(nonlpoly(:, :, ifocal_mod), data.rawdata(:, viewindex));
        end

        % done
        if ~pipeline_onoff
            status.jobdone = true;
        end

    end
end