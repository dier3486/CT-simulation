function [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status)
% recon node, Axial 3D-BP 
% [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status);
% build-in filter and X-upsampling
% cross shot, Z-upsampling
% iteration aligned (filter after Z-alignment)
% call me after reconnode_backprojectionprepare.

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
    [dataflow, prmflow, status] = reconnode_axial3dbackprojectionprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

nodename = status.nodename;
nextnode = status.pipeline.(nodename).nextnode;

% reject while the nextnode is NULL
if pipeline_onoff
    if isempty(nextnode) || strcmpi(nextnode, 'NULL')
        % return;
        status.jobdone = 0;
        status.errorcode = -311;
        status.errormsg = 'The recon node cannot end in last node.';
        return
    end
end

% prio
if pipeline_onoff
    % node prio-step
%     [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    [dataflow, prmflow, status] = BPpriostep(dataflow, prmflow, status);
    % to pass?
    if status.currentjob.topass
        % error or pass
        return;
    end
    % Note: do not passing in last shot
    carrynode = status.currentjob.carrynode;
end

% main
% Nloop
if pipeline_onoff
    % skip looping the shots
    Nloop = 0;
else
    % to loop Nshot+1 times
    Nloop = prmflow.raw.Nshot;
end

% loop the shots
for iloop = 0 : Nloop
    if pipeline_onoff
        [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
            AxialBPKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nodename).data, ...
            dataflow.buffer.(nodename), dataflow.pipepool.(nodename), dataflow.pipepool.(nextnode)(1));
    else
        [dataflow, dataflow.buffer.(nodename)] = AxialBPKernelfuntion(dataflow, dataflow, dataflow.buffer.(nodename));
    end
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done
1;

% Kernel funtion
    function [dataOut, buffer] = AxialBPKernelfuntion(dataOut, dataIn, buffer, currpool, nextpool)
        % The anonymous function is static
        debug = [];

        % parameters
        recon = prmflow.recon;
        Nshot = recon.Nshot;
        Npixel = recon.Npixel;
        Npixel_up = recon.Npixel_up;
        Nslice = recon.Nslice;
        Nviewprot = recon.Nviewprot;
        delta_anglecode = recon.delta_anglecode;
        crosschannel = recon.crosschannel;
        Ncrosschn = recon.Ncrosschn;        % Nchn_up
        Zcrossinterp = recon.Zcrossinterp;
        XY = recon.XY;
        Sxy = recon.activeXY;     % Sxy
        NactiveXY = recon.NactiveXY;  % Nxy
        imagesize = recon.imagesize;
        isXupsampling = recon.upsampling;
        Xupsampgamma = recon.upsampgamma;
        delta_d_up = recon.delta_d_up./recon.SID;
        midchannel_up = recon.midchannel_up;
        Zinterp = recon.Zinterp;
        % Z-upsampling
        Zupsamp = recon.Zupsamp;
        % viewangle
        anglepercode = prmflow.rebin.anglepercode;          % double
        viewangleshift = prmflow.rebin.viewangleshift;      % double
        % iteration
        if isfield(recon, 'iteration_onoff')
            iteration_onoff = recon.iteration_onoff;
        else
            iteration_onoff = false;
        end
        if iteration_onoff
            % to load the iteration parameters
            iteration = prmflow.iteration;
            BV = iteration.BV;
            % noise enhance
            NoiseEnhance = iteration.NoiseEnhance;
        end

        % filters
        basicfilter = recon.filter.basicfilter;
        fftlength = length(basicfilter);
        
        % ini image
        if ~isfield(dataOut, 'image')
            dataOut.image = zeros(imagesize(1)*imagesize(2), 0, 'like', dataIn.rawdata);
        end
        % ini imagehead
        if ~isfield(dataOut, 'imagehead') || isempty(fieldnames(dataOut.imagehead))
            dataOut.imagehead = initialimagehead(0);
        end
        if iteration_onoff && ~isfield(dataOut.imagehead, 'Residual')
            dataOut.imagehead.Residual = zeros(iteration.Niterloop, 0, 'single');
        end

        % ishot
        ishot = buffer.ishot;

        % image number in active
        if ishot == 0 || ishot == Nshot
            if iteration_onoff
                Nimageout = Nslice/2 + 1;
            else
                Nimageout = Nslice/2;
            end
        else
            if iteration_onoff
                Nimageout = Nslice + 2;
            else
                Nimageout = Nslice;
            end
        end
        % image index (in all images)
        if ishot == 0
            imageindex = 1:Nimageout;
        elseif ishot == Nshot
            imageindex = (1:Nimageout) + Nshot*Nslice - Nimageout;
        else
            imageindex = (1:Nimageout) + ishot*Nslice - Nimageout/2;
        end

        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            if ishot < Nshot
                % will get data from dataIn
                NviewIn = plconsol.Index_in(2)-plconsol.Index_in(1)+1;
                index_in = poolindex(currpool, plconsol.Index_in);
                if plconsol.isshotstart
                    startanglecode = double(dataIn.rawhead.Angle_encoder(index_in(1))); % uint32 => int32
                    if ishot > 0
                        % align the viewindex
                        viewshift = double((startanglecode - buffer.startanglecode) / delta_anglecode + Nviewprot/2);
                        % I know buffer.nextshot.poolsize = Nviewprot
                        % to move the WritePoint and ReadPoint
                        buffer.nextshot.ReadStart = mod(buffer.nextshot.ReadStart+viewshift-1, buffer.nextshot.poolsize) + 1;
                        buffer.nextshot.ReadEnd = buffer.nextshot.ReadStart + buffer.nextshot.poolsize - 1;
                        buffer.nextshot.ReadPoint = buffer.nextshot.ReadStart;
                        buffer.nextshot.WritePoint = buffer.nextshot.ReadStart;
                    else % ishot==0
                        buffer.nextshot.ReadEnd = buffer.nextshot.ReadStart + buffer.nextshot.poolsize - 1;
                    end
                    % updata buffer.startanglecode
                    buffer.startanglecode = startanglecode;
                end
                % ishot++ while shot end
                buffer.ishot = buffer.ishot + plconsol.isshotend;
            elseif ishot == Nshot
                if plconsol.isshotstart
                    % start from anywhere is ok
                    buffer.nextshot.ReadStart = 1;
                    buffer.nextshot.ReadEnd = buffer.nextshot.poolsize;
                    buffer.nextshot.ReadPoint = 1;
%                     % close the flag
%                     buffer.lastshotstart = false;
                end
                % only using data in buffer
                NviewIn = min(recon.viewblock, buffer.nextshot.ReadEnd - buffer.nextshot.ReadPoint + 1);
                % to fix the plconsol
                plconsol.isshotend = buffer.nextshot.ReadPoint + NviewIn == buffer.nextshot.ReadEnd + 1;
                % fix the readnumber for last shot (to avoid fake isshotstart)
                status.currentjob.pipeline.readnumber = NviewIn;
                % 
            end
            index_pi = poolindex(buffer.nextshot, [buffer.nextshot.ReadPoint, buffer.nextshot.ReadPoint+NviewIn-1]);
            % to output
            index_out = poolindex(nextpool, plconsol.Index_out);
        else
            startvindex = ishot.*Nviewprot + 1;
            index_in = startvindex:startvindex+Nviewprot-1;
            % align the viewindex
            if ishot < Nshot
                startanglecode = double(dataIn.rawhead.Angle_encoder(index_in(1)));
            else
                startanglecode = double(dataIn.rawhead.Angle_encoder(index_in(1) - Nviewprot));
            end
            if ishot > 0
                prevanglecode = double(dataIn.rawhead.Angle_encoder(index_in(1) - Nviewprot));
                viewshift = double((startanglecode - prevanglecode) / delta_anglecode + Nviewprot/2);
                index_pi = startvindex - Nviewprot + mod((0:Nviewprot-1) + viewshift, Nviewprot);
            end
            NviewIn = Nviewprot;
            % to output
            if iteration_onoff
                index_out = imageindex + ishot*2;
            else
                index_out = imageindex;
            end
            % ishot++
            buffer.ishot = buffer.ishot + 1;
        end

        if recon.couchdirection<0
            % backward couch
            index_leftslice = 1 : (Nslice/2+1)*Npixel;
            index_rightslice = (Nslice/2-1)*Npixel+1 : Nslice*Npixel;
        else
            % forward couch
            index_leftslice = (Nslice/2-1)*Npixel+1 : Nslice*Npixel;
            index_leftslice = reshape(fliplr(reshape(index_leftslice, Npixel, [])), 1, []);
            index_rightslice = 1 : (Nslice/2+1)*Npixel;
            index_rightslice = reshape(fliplr(reshape(index_rightslice, Npixel, [])), 1, []);
        end
        
        % get FBP data       
        switch ishot
            case 0  % first shot
                data0 = dataIn.rawdata(index_leftslice, index_in);
                if isXupsampling
                    data1 = doubleup(reshape(gpuArray(data0), Npixel, []), Xupsampgamma);
                else
                    data1 = reshape(gpuArray(data0), Npixel, []);
                end
                rawdata_ishot = zeros(Ncrosschn, Nslice/2+4, NviewIn, 'like', data1);
                rawdata_ishot(:, 3:end-1, :) = reshape(data1(crosschannel(1):crosschannel(2), :), ...
                    Ncrosschn, Nslice/2+1, NviewIn);
                rawdata_ishot(:, 3:end-1, :) = flipud(rawdata_ishot(:, 3:end-1, :));
                rawdata_ishot(:, 1:2, :) = repmat(rawdata_ishot(:, 3, :), 1, 2);
                rawdata_ishot(:, end, :) = rawdata_ishot(:, end-1, :);
                data_f = zeros(fftlength, Nslice/2+1, NviewIn, 'like', rawdata_ishot);
                % view angle
                AngleEncoder0 = dataIn.rawhead.Angle_encoder(index_in);
                viewangle = cast(double(AngleEncoder0).*anglepercode + viewangleshift - pi/2, 'like', rawdata_ishot);
            case Nshot  % last+1 shot
                if pipeline_onoff
                    data0 = buffer.nextshot.data.rawdata(:, index_pi);
                else
                    data0 = dataIn.rawdata(index_rightslice, index_pi);
                end
                if isXupsampling
                    data1 = doubleup(reshape(gpuArray(data0), Npixel, []), Xupsampgamma);
                else
                    data1 = reshape(gpuArray(data0), Npixel, []);
                end
                rawdata_ishot = zeros(Ncrosschn, Nslice/2+4, NviewIn, 'like', data1);
                rawdata_ishot(:, 2:end-2, :) = reshape(data1(crosschannel(1):crosschannel(2), :), ...
                    Ncrosschn, Nslice/2+1, NviewIn);
                rawdata_ishot(:, 1, :) = rawdata_ishot(:, 2, :);
                rawdata_ishot(:, end-1:end, :) = repmat(rawdata_ishot(:, end-2, :), 1, 2);
                data_f = zeros(fftlength, Nslice/2+1, NviewIn, 'like', rawdata_ishot);
                % view angle
                if pipeline_onoff
                    AngleEncoder0 = buffer.nextshot.data.rawhead.Angle_encoder(index_pi);
                else
                    AngleEncoder0 = dataIn.rawhead.Angle_encoder(index_pi);
                end
                viewangle = cast(double(AngleEncoder0).*anglepercode + viewangleshift + pi/2, 'like', rawdata_ishot);
            otherwise % mid shots
                if pipeline_onoff
                    data0 = [buffer.nextshot.data.rawdata(:, index_pi); dataIn.rawdata(index_leftslice, index_in)];
                else
                    data0 = [dataIn.rawdata(index_rightslice, index_pi); dataIn.rawdata(index_leftslice, index_in)];
                end
                if isXupsampling
                    data1 = doubleup(reshape(gpuArray(data0), Npixel, []), Xupsampgamma);
                else
                    data1 = reshape(gpuArray(data0), Npixel, []);
                end
                rawdata_ishot = zeros(Ncrosschn, Nslice+4, NviewIn, 'like', data1);
                rawdata_ishot(:, 2:end-1, :) = reshape(data1(crosschannel(1):crosschannel(2), :), ...
                    Ncrosschn, Nslice+2, NviewIn);
                rawdata_ishot(:, Nslice/2+3:end-1, :) = flipud(rawdata_ishot(:, Nslice/2+3:end-1, :));
                rawdata_ishot(:, 1, :) = rawdata_ishot(:, 2, :);
                rawdata_ishot(:, end, :) = rawdata_ishot(:, end-1, :);
                data_f = zeros(fftlength, Nslice+2, NviewIn, 'like', rawdata_ishot);
                % view angle
                AngleEncoder0 = dataIn.rawhead.Angle_encoder(index_in);
                viewangle = cast(double(AngleEncoder0).*anglepercode + viewangleshift - pi/2, 'like', rawdata_ishot);
        end

        % save the right slices to private buffer
        if pipeline_onoff
            % I know the nextshot.WritePoint == nextshot.ReadPoint
            if ishot<Nshot
                % copy rawhead
                buffer.nextshot.data = poolhardcopy(buffer.nextshot.data, dataIn, index_pi, index_in, {'rawhead'});
                % copy rawdata of right slices
                buffer.nextshot.data.rawdata(:, index_pi) = dataIn.rawdata(index_rightslice, index_in);
            end
            % move points
            buffer.nextshot.WritePoint = buffer.nextshot.WritePoint + NviewIn;
            buffer.nextshot.ReadPoint = buffer.nextshot.ReadPoint + NviewIn;
        end

        % branch of right part
        if pipeline_onoff && iteration_onoff
            if plconsol.isshotstart
                spaceshift = 0;
            else
                spaceshift = buffer.spaceshift;
            end
            switch ishot
                case 0  % first shot
                    branchpixelindex = Npixel*(Nslice/2+1)+1 : Npixel*(Nslice+2);
                case Nshot  % last+1 shot
                    branchpixelindex = 1 : Npixel*(Nslice/2+1);
                otherwise  % mid shots
                    branchpixelindex = 1 : Npixel*(Nslice+2);
            end
            for isub = 1:iteration.viewsubspace
                % the branch node
                branchnode = plconsol.branchoutput(isub).nextnode;
                branchpoolindex = plconsol.branchoutput(isub).poolindex;
                branchcarrynode = plconsol.branchoutput(isub).carrynode;
                branchcarryindex = plconsol.branchoutput(isub).carryindex;
                % branch subviewshift
                subviewshift = iteration.subviewshift(mod(isub - 1, iteration.viewsubspace) + 1);
                subviewshift = mod(subviewshift + spaceshift - 1, iteration.viewsubspace) + 1;
                % writenum
                writenum = floor((NviewIn - subviewshift) / iteration.viewsubspace) + 1;
                % to return the writenum in branchoutput
                status.currentjob.pipeline.branchoutput(isub).writenumber = writenum;
                status.currentjob.pipeline.branchoutput(isub).newAvail = writenum;
                % copy data0 and AngleEncoder0 to branch pool
                IndexBranch = plconsol.branchoutput(isub).Index;
                IndexBranch(2) = IndexBranch(1) + writenum - 1;
                index_branch = poolindex(dataflow.pipepool.(branchnode)(branchpoolindex), IndexBranch);
                index_subview = subviewshift : iteration.viewsubspace : NviewIn;
                dataflow.pipepool.(branchcarrynode)(branchcarryindex).data.rawdata(branchpixelindex, index_branch) = ...
                    data0(:, index_subview);
                dataflow.pipepool.(branchcarrynode)(branchcarryindex).data.rawhead.Angle_encoder(index_branch) = ...
                    AngleEncoder0(index_subview);

%                 % move the points
%                 [~, dataflow.pipepool.(branchnode)(branchpoolindex)] = ...
%                     movepointsaftercopy([], dataflow.pipepool.(branchnode)(branchpoolindex), [], writenum);
%                 dataflow.pipepool.(branchnode)(branchpoolindex).AvailPoint = ...
%                     dataflow.pipepool.(branchnode)(branchpoolindex).AvailPoint + writenum;
            end
            buffer.spaceshift = mod(spaceshift + (ceil(NviewIn/iteration.viewsubspace)*iteration.viewsubspace - NviewIn), ...
                iteration.viewsubspace);
        elseif iteration_onoff
            % rotate the rawdata
            if ishot>0 && ishot<Nshot
                dataOut.rawdata(:, index_pi+Nviewprot) = dataOut.rawdata(:, index_in);
                for ifield = fieldnames(dataOut.rawhead)'
                    dataOut.rawhead.(ifield{1})(:, index_pi+Nviewprot) = dataOut.rawhead.(ifield{1})(:, index_in);
                end
            end
        end
        
        % prepare for FBP
        switch ishot
            case 0  % first shot
                Zinterp_ishot = Zcrossinterp.t(:, Nslice/2+2:end)-Nslice/2;
                Zinterp_odd_ishot = Zcrossinterp.t_odd(:, Nslice/2+2:end)-Nslice/4;
                Zinterp_even_ishot = Zcrossinterp.t_even(:, Nslice/2+2:end)-Nslice/4;
                Zinterp_gamma_ishot = Zcrossinterp.gamma(:, Nslice/2+2:end);
                chninterp = gpuArray(repmat(single(1:Ncrosschn)', 1, Nslice/2+1));
                interpt_shift = gpuArray(single(-Nslice/2));

                Ztarget = repmat(gpuArray(single(Nslice/2+1 : Nslice/2+Nimageout)), NactiveXY, 1);
                Nslice_shot = gpuArray(Nslice/2);
                ZupMatrix = Zupsamp.ZupMatrix_0;
            case Nshot  % last+1 shot
                Zinterp_ishot = Zcrossinterp.t(:, 1:Nslice/2+1);
                Zinterp_odd_ishot = Zcrossinterp.t_odd(:, 1:Nslice/2+1);
                Zinterp_even_ishot = Zcrossinterp.t_even(:, 1:Nslice/2+1);
                Zinterp_gamma_ishot = Zcrossinterp.gamma(:, 1:Nslice/2+1);
                chninterp = repmat(single(1:Ncrosschn)', 1, Nslice/2+1);
                interpt_shift = gpuArray(single(0));

                Ztarget = repmat(gpuArray(single(Nslice/2-Nimageout+1 : Nslice/2)), NactiveXY, 1);
                Nslice_shot = gpuArray(Nslice/2);
                ZupMatrix = Zupsamp.ZupMatrix_end;
            otherwise
                Zinterp_ishot = Zcrossinterp.t;
                Zinterp_odd_ishot = Zcrossinterp.t_odd;
                Zinterp_even_ishot = Zcrossinterp.t_even;
                Zinterp_gamma_ishot = Zcrossinterp.gamma;
                chninterp = gpuArray(repmat(single(1:Ncrosschn)', 1, Nslice+2));
                interpt_shift = gpuArray(single(0));

                Ztarget = repmat(gpuArray(single(Nslice/2-Nimageout/2+1 : Nslice/2+Nimageout/2)), NactiveXY, 1);
                Nslice_shot = gpuArray(Nslice);
                ZupMatrix = Zupsamp.ZupMatrix_mid;
        end
        
        % cross interps
        for iview = 1:NviewIn
            % 4-point
            data_f(1:Ncrosschn,:, iview) = ...
                interp2(rawdata_ishot(:, 1:2:end, iview), Zinterp_odd_ishot, chninterp, 'linear', 0)./2 + ...
                interp2(rawdata_ishot(:, 2:2:end, iview), Zinterp_even_ishot, chninterp, 'linear', 0)./2 + ...
                interp2(conv2(rawdata_ishot(:,:, iview), [-1 2 -1], 'same'), Zinterp_ishot, chninterp, 'linear', 0).* ...
                Zinterp_gamma_ishot./4;
        end
        % filter
        data_f = reshape(data_f, fftlength, []);
        data_f = ifft(fft(data_f).*basicfilter, 'symmetric');
        data_f = reshape(data_f, fftlength, [], NviewIn);
        % Note: the fftlength could much more than the Npixel_up, in the following calculations only the first Npixel_up values
        % are usable.
        % BP normalize
        BPnormalize = cast(pi/2/Nviewprot, 'like', rawdata_ishot);

        % ini image data to BP
        image_fix = zeros(NactiveXY, Nimageout, 'like', rawdata_ishot);
        for iview = 1:NviewIn
        % for iview = 1: 20 : NviewIn % in test
            % current view angle
            theta = viewangle(iview);

            % (F)BP
            % X-Y to Zeta-Eta
            Eta = -XY(:,1).*sin(theta) + XY(:,2).*cos(theta);
            Zeta = XY(:,1).*cos(theta) + XY(:,2).*sin(theta);

            % to interp on Eta
            Tchn = repmat(Eta./delta_d_up + midchannel_up, 1, Nimageout);

            % Q1
            Tz = interp3(Zinterp.Zeta, Zinterp.Eta, Zinterp.zz, Zinterp.t, ...
                repmat(Zeta, 1, Nimageout), repmat(Eta, 1, Nimageout), Ztarget);
            
            Tz = Tz + interpt_shift;
            Tz = Tz.*(Tz<(Nslice_shot+1) & Tz>0) + (Tz>=(Nslice_shot+1)).*(Nslice_shot+1);
            Tz = Tz.*Zupsamp.Zupsampling + 1;
            D_UP = data_f(:, :, iview)*ZupMatrix;

            image_fix = image_fix + interp2(D_UP, Tz, Tchn, 'linear', 0);
        end
        % normalize by pi/2/Nviewprot
        image_fix = image_fix.*BPnormalize;
        % add to output
        if size(dataOut.image, 2) < max(index_out)
            dataOut.image(:, max(index_out)) = 0;
        end
        dataOut.image(Sxy(:), index_out) = dataOut.image(Sxy(:), index_out) + image_fix;

        % image head and TV
        if ~pipeline_onoff || plconsol.isshotend
            dataOut.imagehead = reconimagehead(dataOut.imagehead, recon, imageindex, index_out, ishot);
            % TV
            if iteration_onoff
                % get images to do TV
                u = cast(reshape(dataOut.image(:, index_out), imagesize(2), imagesize(1), Nimageout), 'like', rawdata_ishot);
                % Bregman TV
                u = u.*(1 + 1i);
                [u, b] = BregmanTV3D_axialiter(u, BV.mu_BV0, BV.Cl, [], BV.Crange, BV.maxNiter, BV.tol2, ...
                    BV.Zbalance);
                mu_BV = BV.mu_BV0./(abs(real(u) - imag(u)).*BV.mu_adap + 1);
                u = BregmanTV3D_axialiter(u, mu_BV, BV.Cl, b, BV.Crange, BV.maxNiter, BV.tol, BV.Zbalance);
                % noise enhance
                u = enhancemix_iter(u, NoiseEnhance.wlevel, NoiseEnhance.mixwidth, NoiseEnhance.mixalpha, NoiseEnhance.mixbeta);
                % or, u = real(u) + (real(u)-imag(u)).*(-1/2+1/2i);
                % output
                u = reshape(u, [], Nimageout);
                dataOut.image(:, index_out) = gather(u);
            end
        end

        % jobdone fix
        if pipeline_onoff
            if ishot==Nshot-1 && status.jobdone==1
                % the Nshot-1 shot is special due to need to wake up the last shot.
                status.jobdone = 2;
                % keep waking to do the last shot
            elseif ishot==Nshot
                % for the last shot the plc
                if plconsol.isshotend
                    % job done
                    status.currentjob.pipeline.writenumber = Nimageout;
                    status.currentjob.pipeline.newAvail = Nimageout;
                    status.jobdone = 1;
                else
                    % call me gain
                    status.jobdone = 7;        
                end
            end
        else
            status.jobdone = true;
        end
        % BP Kernelfuntion END
    end

end


function [dataflow, prmflow, status] = BPpriostep(dataflow, prmflow, status)
% Axial

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;
% I know it was rejected while the nextnode is NULL.
pipeprm = nodeprm.pipeline;

Nshot = prmflow.recon.Nshot;
Nslice = prmflow.recon.Nslice;
Nviewprot = prmflow.recon.Nviewprot;
Nimage = prmflow.recon.Nimage;
% iteration
if isfield(prmflow.recon, 'iteration_onoff')
    iteration_onoff = prmflow.recon.iteration_onoff;
else
    iteration_onoff = false;
end

ishot = dataflow.buffer.(nodename).ishot;
if ishot==0 || ishot==Nshot
    if ~iteration_onoff
        Nimage_ishot = Nslice/2;
    else
        Nimage_ishot = Nslice/2 + 1;
    end
else
    if ~iteration_onoff
        Nimage_ishot = Nslice;
    else
        Nimage_ishot = Nslice + 2;
    end
end

% ini the flag to run the poststep
status.currentjob.torunpoststep = true;

% prio step #2, the next node and carry node
status.currentjob.nextnode = nextnode;
if ~isempty(dataflow.pipepool.(nextnode))
    if dataflow.pipepool.(nextnode)(1).iscarried
        carrynode = dataflow.pipepool.(nextnode)(1).carrynode;
    else
        carrynode = nextnode;
    end
else
    carrynode = nextnode;
end
status.currentjob.carrynode = carrynode;

% prio step #1, the prio steps of the pipepools
% pass due to nextpool is stucked
nextWritableNumber = dataflow.pipepool.(nextnode)(1).poolsize - dataflow.pipepool.(nextnode)(1).WritePoint + 1;
if (nextWritableNumber < Nimage_ishot) || dataflow.pipepool.(nextnode)(1).WriteStuck
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 6;
    % keep waking
    % to pass the kernel function
    status.currentjob.topass = true;
    return;
end

% I know the nodetype is A-H.1.G or H-H.1.G

% check shot start 1
% status.currentjob.pipeline.isshotstart = ...
%     dataflow.pipepool.(nodename)(1).ReadPoint == dataflow.pipepool.(nodename)(1).ReadStart;
status.currentjob.pipeline.isshotstart = dataflow.pipepool.(nextnode)(1).isshotstart;
% close it
dataflow.pipepool.(nextnode)(1).isshotstart = false;

% ini next pool
if status.currentjob.pipeline.isshotstart
        dataflow.pipepool.(nextnode)(1).ReadStart = 1;
        dataflow.pipepool.(nextnode)(1).ReadEnd = Nimage_ishot;
        dataflow.pipepool.(nextnode)(1).WriteStart = 1;
        dataflow.pipepool.(nextnode)(1).WriteEnd = Nimage_ishot;
        dataflow.pipepool.(nextnode)(1).ReadPoint = 1;
        dataflow.pipepool.(nextnode)(1).WritePoint = 1;
        dataflow.pipepool.(nextnode)(1).AvailPoint = 0;

    % shot start in iteration recon
    dataflow.pipepool.(nextnode)(1).ReadStart = 1;
    dataflow.pipepool.(nextnode)(1).ReadEnd = Nimage_ishot;
    dataflow.pipepool.(nextnode)(1).WriteStart = 1;
    dataflow.pipepool.(nextnode)(1).WriteEnd = Nimage_ishot;
    dataflow.pipepool.(nextnode)(1).ReadPoint = 1;
    dataflow.pipepool.(nextnode)(1).WritePoint = 1;
    dataflow.pipepool.(nextnode)(1).AvailPoint = 0;
    if iteration_onoff
        % branchs
        viewsubspace = prmflow.recon.viewsubspace;
        for ii = 1: viewsubspace
            % the branch node
            branchnode = status.pipeline.(nodename).branchnextnodes{ii};
            branchpoolindex = status.pipeline.(nodename).branchnextpoolindex(ii);
            if dataflow.pipepool.(branchnode)(branchpoolindex).WriteStuck
                % stuck due to the branchpool
                status.currentjob.pipeline.readnumber = 0;
                status.currentjob.pipeline.writenumber = 0;
                status.currentjob.pipeline.newAvail = 0;
                status.jobdone = 6;
                % to pass the kernel function
                status.currentjob.topass = true;
                % withdraw the initial of next pool
                dataflow.pipepool.(nextnode)(1).isshotstart = true;
                return;
            end
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadEnd = Nviewprot / viewsubspace;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteEnd = Nviewprot / viewsubspace;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadPoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WritePoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).AvailPoint = 0;
        end
    end
end

% trace
if status.currentjob.pipeline.isshotstart && status.debug.pooltrace_onoff
    if ~iteration_onoff
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'priostep');
    else
        dataflow.pipepool = pooltrace(dataflow.pipepool, nodename, nextnode, 'priostep', ...
            status.pipeline.(nodename).branchnextnodes, status.pipeline.(nodename).branchnextpoolindex);
    end
end

% check shot end 1 (the input data has reached the end)
Ishotend1 = dataflow.pipepool.(nodename)(1).WritePoint == dataflow.pipepool.(nodename)(1).WriteEnd + 1;

% minlimit
status.currentjob.pipeline.minlimit = pipeprm.inputminlimit;

% n (currAvailNumber), m
currAvailNumber = dataflow.pipepool.(nodename)(1).AvailPoint - dataflow.pipepool.(nodename)(1).ReadPoint + 1;
n = min(pipeprm.inputmaxlimit, currAvailNumber);
m = Nimage_ishot;

if ishot~=Nshot && (n < status.currentjob.pipeline.minlimit && ~Ishotend1)
    % not enough input views
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 3;
    % to pass the kernel function
    status.currentjob.topass = true;
    % withdraw the initial of next pool
    dataflow.pipepool.(nextnode)(1).isshotstart = status.currentjob.pipeline.isshotstart;
    return;
end
% we skipped to check if the space of the branchpools are enough which are believed in that.

% check shot end 2 (will output the shot end)
if ishot < Nshot
    status.currentjob.pipeline.isshotend = Ishotend1 && n == currAvailNumber;
else
    status.currentjob.pipeline.isshotend = false;
    % while the ishot==Nshot plz check the shotend in kernal function
end

% Index_in, Index_out
status.currentjob.pipeline.Index_in = ...
    [dataflow.pipepool.(nodename)(1).ReadPoint dataflow.pipepool.(nodename)(1).ReadPoint+n-1];
status.currentjob.pipeline.Index_out = ...
    [dataflow.pipepool.(nextnode)(1).WritePoint dataflow.pipepool.(nextnode)(1).WritePoint+m-1];
% Do2i (not use)
status.currentjob.pipeline.Do2i = 0;
% Nexpand
status.currentjob.pipeline.Nexpand = 0;

% readnumber
status.currentjob.pipeline.readnumber = n;

% writenumber, newAvail
if status.currentjob.pipeline.isshotend
    status.currentjob.pipeline.writenumber = Nimage_ishot;
    status.currentjob.pipeline.newAvail = Nimage_ishot;
else
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
end

% Done
if status.currentjob.pipeline.isshotend
    status.jobdone = 1;
elseif n < currAvailNumber
    status.jobdone = 7;
else
    status.jobdone = 4;
end

% if n < currAvailNumber
%     status.jobdone = 7;
% else
%     status.jobdone = 1;
% end

% branchs
if iteration_onoff
    viewsubspace = prmflow.iteration.viewsubspace;
    status.currentjob.pipeline.branchoutput(viewsubspace) = struct();
    for isub = 1 : viewsubspace
        % the branch node
        branchnode = status.pipeline.(nodename).branchnextnodes{isub};
        branchpoolindex = status.pipeline.(nodename).branchnextpoolindex(isub);
        if dataflow.pipepool.(branchnode)(branchpoolindex).iscarried
            carrynode = dataflow.pipepool.(branchnode)(branchpoolindex).carrynode;
            carryindex = dataflow.pipepool.(branchnode)(branchpoolindex).carryindex;
        else
            carrynode = branchnode;
            carryindex = branchpoolindex;
        end
        status.currentjob.pipeline.branchoutput(isub).nextnode = branchnode;
        status.currentjob.pipeline.branchoutput(isub).poolindex = branchpoolindex;
        status.currentjob.pipeline.branchoutput(isub).carrynode = carrynode;
        status.currentjob.pipeline.branchoutput(isub).carryindex = carryindex;
        status.currentjob.pipeline.branchoutput(isub).writenumber = nan;    % to be set
        status.currentjob.pipeline.branchoutput(isub).newAvail = nan;       % to be set
        status.currentjob.pipeline.branchoutput(isub).Index = ...
            [dataflow.pipepool.(branchnode)(branchpoolindex).WritePoint -1];   % to be set, or skip
    end
end



end