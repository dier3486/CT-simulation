function [dataflow, prmflow, status] = reconnode_Antiwindmill(dataflow, prmflow, status)
% recon node, anti windmill artifact on image
% TV based method
% [dataflow, prmflow, status] = reconnode_Antiwindmill(dataflow, prmflow, status);
% suggusted position: unfixed.

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
    [dataflow, prmflow, status] = reconnode_antiwindmillprepare(dataflow, prmflow, status);
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
    [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
        AntiwindmillKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nextnode)(1), ...
        dataflow.buffer.(nodename));
else
    dataflow = AntiwindmillKernelfuntion(dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function [dataOut, buffer] = AntiwindmillKernelfuntion(dataOut, nextpool, buffer)
        % The anonymous function is static
        debug = [];

        % no images?
        if ~isfield(dataOut, 'image')
            return;
        end
        
        % corr table
        awindcorr = prmflow.corrtable.(nodename);
        TVmu = awindcorr.TVmu;
        TVlambda = awindcorr.TVlambda; 
        TVlogC = awindcorr.TVlogC;
        TVCrange = awindcorr.TVCrange;
        TVNiter = awindcorr.TVNiter; 
        TVtol = awindcorr.TVtol;
        fixsigma = awindcorr.fixsigma;
        fixlimit = awindcorr.fixlimit;

        nodeprm = prmflow.pipe.(nodename);

        imagesize = prmflow.image.imagesize;
        % GPU on/off
        GPUonoff = status.currentjob.GPUdevice > 0;

        % is real?
        flag_real = isreal(dataOut.image);
%         if flag_real
%             nodeprm = everything2real(nodeprm);
%         end

        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            imagerely = nodeprm.pipeline.viewrely;
            imagerely_in = [imagerely(2) imagerely(2)];
            if plconsol.isshotstart
                imagerely_in(1) = 0;
            end
            if plconsol.isshotend
                imagerely_in(2) = 0;
            end
            index_in = poolindex(nextpool, plconsol.Index_out + imagerely_in);
            % NimageIn = length(index_in);
            index_out = poolindex(nextpool, plconsol.Index_out);
            NimageOut = length(index_out);
            index_pick = 1:NimageOut;
        else
            imagerely = [0 0];
            NimageOut = prmflow.recon.Nimage;
            index_in = 1 : NimageOut;
            index_out = index_in;
            index_pick = index_in;
        end

        % get data
        image0 = reshape(dataOut.image(:, index_in), imagesize(2), imagesize(1), []);

        % Gauss blur
        if abs(awindcorr.Gblur) > eps
            image0 = gaussblur(image0, awindcorr.Gblur);
            % I know the gaussblur can treat complex data.
        end
        
        % down sampling
        if awindcorr.downsample > 1
            image0 = imagedown(image0, awindcorr.downsample);
        end

        if pipeline_onoff
            if ~plconsol.isshotstart
                image0 = cat(3, buffer.imagebound, image0);
                Nbnd = size(buffer.imagebound, 3);
                index_pick = index_pick + Nbnd - imagerely(2);
            end
            if ~plconsol.isshotend
                buffer.imagebound = image0(:,:, max(end - sum(imagerely) + 1, 1) : end);
            end
        end

        % expanded TV
        DIM = 3;
        if flag_real
            imageTV = GlogTV1D(image0, TVmu, TVlambda, TVlogC, [], DIM, TVCrange, TVNiter, TVtol);
        else
            imageTV = complex(GlogTV1D(real(image0), real(TVmu), real(TVlambda), real(TVlogC), [], DIM, real(TVCrange), ...
                real(TVNiter), real(TVtol)), ...
                GlogTV1D(imag(image0), imag(TVmu), imag(TVlambda), imag(TVlogC), [], DIM, imag(TVCrange), ...
                imag(TVNiter), imag(TVtol)) );
        end

        % diff
        image0 = imageTV(:,:, index_pick) - image0(:,:, index_pick);

        % limit
        if flag_real
            image0 = hardentanh(image0, fixsigma, fixlimit);
        else
            image0 = complex(hardentanh(real(image0), real(fixsigma), real(fixlimit)), ...
                hardentanh(imag(image0), imag(fixsigma), imag(fixlimit)) );
        end

        % up sampling
        if awindcorr.downsample > 1
            image0 = imageupnGPU(image0, awindcorr.downsample);
            image0 = image0(1:imagesize(2), 1:imagesize(1), :);
        end
        
        % reshape
        image0 = reshape(image0, [], NimageOut);

%         % for sparse image
%         if isfield(dataOut, 'Sreconact')  % sparse image, we shall have an object type flag in prmflow.image
%             image0 = image0.*dataOut.Sreconact(:, index_out);
%         end

        % fix
        if pipeline_onoff
            dataOut.image(:, index_out) = dataOut.image(:, index_out) + image0;
        else
            dataOut.image(:, index_out) = dataOut.image(:, index_out) + gather(image0);
        end

        if ~pipeline_onoff
            status.jobdone = true;
        end
        % AntiwindmillKernelfuntion END
    end

end

