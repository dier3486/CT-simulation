function [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status)
% recon node, interp for bad channels
% [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_badchannelprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

nodename = status.nodename;
nextnode = status.pipeline.(nodename).nextnode;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
    carrynode = status.currentjob.carrynode;
end

% main
if pipeline_onoff
    dataflow.pipepool.(carrynode)(1).data = badchannelKernelfuntion(dataflow.pipepool.(carrynode)(1).data, ...
        dataflow.pipepool.(nextnode)(1));
else
    dataflow = badchannelKernelfuntion(dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function dataOut = badchannelKernelfuntion(dataOut, nextpool)
        % The anonymous function is static
        debug = [];
        
        % parameters
        nodeprm = prmflow.pipe.(nodename);
        badindex = nodeprm.badindex;
        Nbadchennel = nodeprm.Nbadchennel;
        interpindex = nodeprm.interpindex;
        interpalpha = nodeprm.interpalpha;
        
        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(nextpool, plconsol.Index_out);
        else
            index_out = 1 : prmflow.raw.Nview;
        end

        % we nolonger to check the NaN

        % loop the bad channels
        for ii = 1:Nbadchennel
            % interp
            dataOut.rawdata(badindex(ii), index_out) = interpalpha(ii, :) * dataOut.rawdata(interpindex(ii, :), index_out);
        end
        
        % job done
        if ~pipeline_onoff
            status.jobdone = true;
        end
    end
end