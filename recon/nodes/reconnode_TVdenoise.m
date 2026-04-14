function [dataflow, prmflow, status] = reconnode_TVdenoise(dataflow, prmflow, status)
% post recon node, TVdenoise (after Backprojection)
% [dataflow, prmflow, status] = reconnode_TVdenoise(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_tvdenoiseprepare(dataflow, prmflow, status);
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
    % CUDA on/off
    CUDAonoff = prmflow.pipe.(nodename).pipeline.CUDAonoff;
end

% main
if pipeline_onoff
    if CUDAonoff
        [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
            DenoiseCUDAfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename), ...
            dataflow.pipepool.(nextnode)(1), status.currentjob.pipeline);
    else
        dataflow.pipepool.(carrynode)(1).data = ...
            DenoiseKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nextnode)(1));
    end
else
    dataflow = DenoiseKernelfuntion(dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function dataOut = DenoiseKernelfuntion(dataOut, nextpool)
        % The anonymous function is static
        debug = [];

        % no images?
        if ~isfield(dataOut, 'image')
            return;
        end

        % prm
        pipeprm = prmflow.pipe.(nodename).pipeline;
        imageprm = prmflow.image;
        imagesize = prmflow.recon.imagesize;


        % GPU on/off
        GPUonoff = status.currentjob.GPUdevice > 0;

        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            isshotstart = plconsol.isshotstart;
            % index_in = poolindex(nextpool, plconsol.Index_in);
            index_out = poolindex(nextpool, plconsol.Index_out);
            IndexIn = plconsol.Index_in;
            if ~isshotstart
                IndexIn(1) = plconsol.Index_out(1) - pipeprm.viewrely(1);
            end
            index_in = poolindex(nextpool, IndexIn);
            
            NimageIn = length(index_in);
            NimageOut = length(index_out);
        else
            isshotstart = true;
            NimageOut = prmflow.recon.Nimage;
            index_out = 1 : NimageOut;
            NimageIn = NimageOut;
            index_in = index_out;
        end

        % get data
        imgfix = real(reshape(dataOut.image(:, index_in), imagesize(2), imagesize(1), NimageIn));
        % only fix the real part of the images

        % fill nan (can skip)
        imgfix(isnan(imgfix)) = 0;
        
        % TVopts
        TVopts = imageprm.TVopts;
        if strcmp(TVopts.DIM, '3D') && ~isshotstart && pipeprm.viewrely(1)>0
            TVopts.LeftBoundary = true;
        end

        % BegmanTV
        TVoptsC = namedargs2cell(TVopts);
        imgfix = BregmanTV(imgfix, imageprm.TV_mu, imageprm.TV_Cl, [], [], TVoptsC{:});

        % out
        imgfix = reshape(imgfix, [], NimageIn);
        index_fixout = (1 : NimageOut) + index_out(1) - index_in(1);
        dataOut.image(:, index_out) = imgfix(:, index_fixout) + imag(dataOut.image(:, index_out)).*1i;

    end
end