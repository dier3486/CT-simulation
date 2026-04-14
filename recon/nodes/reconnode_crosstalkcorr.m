function [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status)
% recon node, crosstalk correction
% [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status);
% suggested position: before Beamharden/Nonliear, after Reference (and other correctons need to preset to Crosstalk).

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
    [dataflow, prmflow, status] = reconnode_crosstalkprepare(dataflow, prmflow, status);
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
if pipeline_onoff
    nextnode = status.currentjob.nextnode;
    carrynode = status.currentjob.carrynode;
    dataflow.pipepool.(carrynode).data = ...
        crosstalkcorrKernelfuntion(dataflow.pipepool.(nextnode), dataflow.pipepool.(carrynode).data);
else
    dataflow = crosstalkcorrKernelfuntion([], dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

    % Kernel funtion
    function data = crosstalkcorrKernelfuntion(nextpool, data)
        % The anonymous function is static
        debug = [];

        % parameters to use in prmflow
        Npixel = prmflow.raw.Npixel;
        Nslice = prmflow.raw.Nslice;

        % parameters set in pipe
        crossprm = prmflow.pipe.(status.nodename);
        weight = crossprm.weight;
        istointensity = crossprm.istointensity;

        % calibration table
        crscorr = prmflow.corrtable.(status.nodename);
        crsorder = crscorr.order;
        % DFS
        Nfocal = prmflow.raw.Nfocal;
        crsval_org = reshape(crscorr.main, [], crsorder, crscorr.focalnumber);

        % GPU on/off
        GPUonoff = status.currentjob.GPUdevice > 0;

        % data classes
        if GPUonoff
            dataclass_single = ones(1, 'single', 'gpuArray');
            dataclass_raw = gpuArray(ones(1, 'like', data.rawdata));
        else
            dataclass_single = ones(1, 'single');
            dataclass_raw = ones(1, 'like', data.rawdata);
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

        % step1, to intensity
        if istointensity
            data.rawdata(:, index_out) = 2.^(-data.rawdata(:, index_out));
        end

        % step2, correct
        for ifocal = 1:Nfocal
            viewindex = index_out(ifocal : Nfocal : Nview);
            ifocal_mod = mod((ifocal-1), crscorr.focalnumber) + 1;
            crsval = crsval_org(:, :, ifocal_mod);
            switch crsorder
                case 1
                    % odd-symmetric style
                    % the correction operator is a tridiagonal matrix [-crsval;  1; crsval];
                    % rawfix
                    rawfix = zeros(Npixel*Nslice, length(viewindex), 'like', dataclass_raw);
                    rawfix(2:end-1, :) = (data.rawdata(3:end, viewindex) - data.rawdata(1:end-2, viewindex)).*crsval(2:end-1);
                    % add to rawdata
                    data.rawdata(:, viewindex) = data.rawdata(:, viewindex) + rawfix.*weight;
                case 3
                    % 0
                    rawfix = data.rawdata(:, viewindex).*crsval(:,2);
                    % -1
                    rawfix(2:end, :) = rawfix(2:end, :) + data.rawdata(1:end-1, viewindex).*crsval(2:end, 1);
                    % 1
                    rawfix(1:end-1, :) = rawfix(1:end-1, :) + data.rawdata(2:end, viewindex).*crsval(1:end-1, 3);
                    % add to rawdata
                    data.rawdata(:, viewindex) = data.rawdata(:, viewindex) + rawfix.*weight;
                case 5
                    % 0
                    rawfix = data.rawdata(:, viewindex).*crsval(:,3);
                    % -2
                    rawfix(3:end, :) = rawfix(3:end, :) + data.rawdata(1:end-2, viewindex).*crsval(3:end, 1);
                    % -1
                    rawfix(2:end, :) = rawfix(2:end, :) + data.rawdata(1:end-1, viewindex).*crsval(2:end, 2);
                    % 1
                    rawfix(1:end-1, :) = rawfix(1:end-1, :) + data.rawdata(2:end, viewindex).*crsval(1:end-1, 4);
                    % 2
                    rawfix(1:end-2, :) = rawfix(1:end-2, :) + data.rawdata(3:end, viewindex).*crsval(1:end-2, 5);
                    % add to rawdata
                    data.rawdata(:, viewindex) = data.rawdata(:, viewindex) + rawfix.*weight;
                otherwise
                    error('Crosstalk correction not support order %d.', crsorder);
            end
        end

        if istointensity
            % min cut
            minval = cast(2^-32, 'like', dataclass_single);
            % log2
            data.rawdata(:, index_out) = log2opterators(data.rawdata(:, index_out), minval);
        end

        % done
        if ~pipeline_onoff
            status.jobdone = true;
        end

    end
end

function rawdata = log2opterators(rawdata, minval)

rawdata(rawdata < minval) = minval;
% log2
rawdata = -log2(rawdata);

end