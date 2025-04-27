function [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status)
% recon node, filter design and conv
% [dataflow, prmflow, status] = reconnode_filter(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_filterprepare(dataflow, prmflow, status);
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
    dataflow.pipepool.(carrynode)(1).data = filterKernelfuntion(dataflow.pipepool.(carrynode)(1).data, ...
        dataflow.pipepool.(nodename)(1).data,...
        dataflow.pipepool.(nextnode)(1), dataflow.pipepool.(nodename)(1));
else
    dataflow = filterKernelfuntion(dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function dataOut = filterKernelfuntion(dataOut, dataIn, nextpool, currpool)
        % The anonymous function is static
        debug = [];

        % paramters
        nodeprm = prmflow.pipe.(nodename);
        isXupsampling = nodeprm.upsampling;
        Nslice = prmflow.recon.Nslice;
        Npixel = prmflow.recon.Npixel;
        if isXupsampling
            Npixel_up = prmflow.recon.Npixel_up;
            mid_u = prmflow.recon.midchannel_up;
        else
            Npixel_up = Npixel;
            mid_u = prmflow.recon.midchannel;
        end

        % filter
        basicfilter = prmflow.recon.filter.basicfilter;
        Hlen = prmflow.recon.filter.Hlen;

        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_out = poolindex(nextpool, plconsol.Index_out);
            NviewOut = length(index_out);
            kernellevel = nodeprm.pipeline.kernellevel;
            if kernellevel ~= 0
                index_in = poolindex(currpool, plconsol.Index_in);
            end
        else
            index_out = 1 : prmflow.rebin.NviewOut;
            NviewOut = prmflow.rebin.NviewOut;
            kernellevel = 0;
            % index_in = index_out;
        end
        
        % get data
        A = zeros(Hlen, Nslice*NviewOut, 'like', dataOut.rawdata);
        if kernellevel == 0
            if isXupsampling % I know, pipeline_onoff==false
                A(1:2:Npixel_up, :) = reshape(dataOut.rawdata(:, index_out), Npixel, []);
            else
                A(1:Npixel, :) = reshape(dataOut.rawdata(:, index_out), Npixel, []);
            end
        else % kernellevel==1, mostly isXupsampling==true
            if isXupsampling
                A(1:2:Npixel_up, :) = reshape(dataIn.rawdata(:, index_in), Npixel, []);
            else
                A(1:Npixel, :) = reshape(dataIn.rawdata(:, index_in), Npixel, []);
            end
        end

        % upsampling
        if isXupsampling
            A(1:Npixel_up, :) = doubleup2(A(1:Npixel_up, :), nodeprm.upsampgamma);
        end
        
        % fill up
        if isfield(nodeprm, 'fillup') && nodeprm.fillup
            % to fill up the data for off-ISOcenter
            if kernellevel == 0
                if isfield(dataOut.rawhead, 'refblock')
                    blkvindex = any(dataOut.rawhead.refblock(:, index_out), 1);
                else
                    blkvindex = [];
                end
            else
                if isfield(dataIn.rawhead, 'refblock')
                    blkvindex = any(dataIn.rawhead.refblock(:, index_in), 1);
                else
                    blkvindex = [];
                end
            end
            A = reshape(A, Hlen, Nslice, []);
            for ii = 1:Nslice
                [A(:, ii, :), n_left] = translatefillup(squeeze(A(1:Npixel_up, ii, :)), Hlen, mid_u, blkvindex);
            end
            A = reshape(A, Hlen, Nslice*NviewOut);
            % the fill up is only used to observe the off-ISOcenter water phantom in calibration
        else
            n_left = 0;
        end

        % isreal
        isreal_flag = isreal(dataOut.rawdata);

        % conv
        % fft
        A = fft(A);
        % timesymmetric
        A = A.*basicfilter;
        % ifft
        if isreal_flag
            A = ifft(A, 'symmetric');
        else
            A = ifft(A);
        end
        
        % kick filled zero
        if pipeline_onoff
            dataOut.rawdata(:, index_out) = reshape(A((1:Npixel_up)+n_left, :), Npixel_up*Nslice, []);
        else
            dataOut.rawdata = reshape(A((1:Npixel_up)+n_left, :), Npixel_up*Nslice, []);
        end

        % done
        if ~pipeline_onoff
            status.jobdone = true;
        end
    end

end