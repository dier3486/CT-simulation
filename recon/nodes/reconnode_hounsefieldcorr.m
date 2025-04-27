function [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status)
% recon node, Hounsefield Units correction
% [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status);

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

% no prepare of Hounsefield correction

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.jobdone==0 || status.jobdone>=3
        % error or pass
        return;
    end
end

% main
if pipeline_onoff
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
    dataflow.pipepool.(carrynode).data = ...
        hounsefieldcorrKernelfuntion(dataflow.pipepool.(nextnode), dataflow.pipepool.(carrynode).data);
else
    dataflow = hounsefieldcorrKernelfuntion([], dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function data = hounsefieldcorrKernelfuntion(nextpool, data)
        % parameters to use in prmflow
        Npixel = prmflow.raw.Npixel;
        Nslice = prmflow.raw.Nslice;
        % parameters set in pipe
        HCprm = prmflow.pipe.(status.nodename);
        if isfield(HCprm, 'HCscale')
            HCscale = HCprm.HCscale;
        else
            HCscale = 1000;
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

        % scale 1000
        data.rawdata(:, index_out) = data.rawdata(:, index_out).*HCscale;

        % calibration table
        if isfield(prmflow.corrtable, status.nodename)
            % the Hounsefield calibration is not always necessary.
            HUcorr = prmflow.corrtable.(status.nodename);
            data.rawdata(:, index_out) = reshape(reshape(data.rawdata(:, index_out), Npixel, Nslice, Nview).*HUcorr.main(:)'...
                , [], Nview);
        end

        % done
        if ~pipeline_onoff
            status.jobdone = true;
        end
    end
end