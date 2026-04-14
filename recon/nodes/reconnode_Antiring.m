function [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status)
% recon node, anti ring on image
% [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);
% suggusted position: after Antiwindmill (should be better), and before other recon-postprocess nodes.

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
    [dataflow, prmflow, status] = reconnode_antiringprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

nodename = status.nodename;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
    nextnode = status.pipeline.(nodename).nextnode;
    carrynode = status.currentjob.carrynode;
end

% main
if pipeline_onoff
    dataflow.pipepool.(carrynode)(1).data = ...
        AntiringKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nextnode)(1));
else
    dataflow = AntiringKernelfuntion(dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function dataOut = AntiringKernelfuntion(dataOut, nextpool)
        % The anonymous function is static
        debug = [];
        
        % no images?
        if ~isfield(dataOut, 'image')
            return;
        end

        % corr table
        aringcorr = prmflow.corrtable.(nodename);
        alpha = aringcorr.alpha;
        Ntheta = aringcorr.Ntheta;
        Crange = aringcorr.Crange;
        rtheta_evenodd = aringcorr.rtheta_evenodd;
        ringcut = aringcorr.ringcut;
        ringmover = aringcorr.ringmover;
        ringsection = aringcorr.ringsection;
        % parameters
        imagesize = prmflow.image.imagesize;
        voxelsize = prmflow.image.voxelsize;
        d_radius = prmflow.image.delta_radius;
        Rfov = prmflow.image.effFOV/2;
        % GPU on/off
        GPUonoff = status.currentjob.GPUdevice > 0;

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(nextpool, plconsol.Index_out);
            NimageOut = length(index_out);
        else
            NimageOut = prmflow.recon.Nimage;
            index_out = 1 : NimageOut;
        end
        
        % get data
        imgfix = real(reshape(dataOut.image(:, index_out), imagesize(2), imagesize(1), []));
        % only fix the real part of the images
        imagecenter = dataOut.imagehead.imagecenter(:, index_out);

        % fill nan (can skip)
        % imgfix = fillmissing(imgfix, 'constant', 0);  % the fillmissing is very slow
        imgfix(isnan(imgfix)) = 0;

        if GPUonoff
            % putinGPU
            [imagecenter, voxelsize, d_radius] = ...
                putinGPU(imagecenter, voxelsize, d_radius);
        end
        centerfix = imagecenter(1:2, :)'./voxelsize;
        
        % anti ring (only on real part of the images)
        imgfix = antiringonimage(imgfix, centerfix, Crange, ...
            Ntheta, d_radius, Rfov, rtheta_evenodd, ringcut, ringmover, ringsection);

        % reshape and  *alpha
        imgfix = reshape(imgfix, [], NimageOut).*alpha;

%         % for sparse image
%         if isfield(dataOut, 'Sreconact')  % sparse image, we shall have an object type flag in prmflow.image
%             imgfix = imgfix.*dataOut.Sreconact(:, index_out);
%         end

        % add to image
        dataOut.image(:, index_out) = dataOut.image(:, index_out) - imgfix;

        if ~pipeline_onoff
            status.jobdone = true;
        end
        % AntiringKernelfuntion END
    end

end