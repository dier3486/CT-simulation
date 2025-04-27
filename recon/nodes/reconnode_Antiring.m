function [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status)
% recon node, anti ring on image
% [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);

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

        % GPU?
        if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
            GPUonoff = true;
        else
            GPUonoff = false;
        end
        
        % no images?
        if ~isfield(dataOut, 'image')
            return;
        end

        % parameters
        nodeprm = prmflow.pipe.(nodename);
        alpha = nodeprm.alpha;
        Ntheta = nodeprm.Ntheta;
        Crange = nodeprm.Crange;
        Cnorm = nodeprm.Cnorm;
        rtheta_evenodd = nodeprm.rtheta_evenodd;
        restcutoff = nodeprm.restcutoff;
        ringfilter = nodeprm.ringfilter;
        ringsection = nodeprm.ringsection;
        sectmethod = nodeprm.sectmethod;
        
        imagesize = prmflow.recon.imagesize;
        voxelsize = prmflow.recon.voxelsize;
        d_radius = prmflow.recon.delta_d/voxelsize;
                
%         % is real?
%         flag_real = isreal(dataOut.image);

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(nextpool, plconsol.Index_out);
            NimageOut = length(index_out);
        else
            NimageOut = prmflow.recon.Nimage;
            index_out = 1 : NimageOut;
        end
        
        % get data
        if GPUonoff
            imgfix = gpuArray(reshape(dataOut.image(:, index_out), imagesize(2), imagesize(1), []));
            imagecenter = gpuArray(dataOut.imagehead.imagecenter(:, index_out));
            % putinGPU
            [voxelsize, Crange, Cnorm, Ntheta, d_radius, restcutoff, ringfilter] = ...
                putinGPU(voxelsize, Crange, Cnorm, Ntheta, d_radius, restcutoff, ringfilter);
        else
            imgfix = reshape(dataOut.image(:, index_out), imagesize(2), imagesize(1), []);
            imagecenter = dataOut.imagehead.imagecenter(:, index_out);
        end
        centerfix = imagecenter(1:2, :)'./voxelsize;
        
        % anti ring (only on real part of the images
        imgfix = antiringonimage(real(imgfix), centerfix, Crange, Cnorm, Ntheta, d_radius, rtheta_evenodd, restcutoff, ...
            ringfilter, ringsection, sectmethod);

        % add to image
        dataOut.image = dataOut.image - reshape(gather(imgfix), [], NimageOut).*alpha;

        if ~pipeline_onoff
            status.jobdone = true;
        end
        % AntiringKernelfuntion END
    end

end