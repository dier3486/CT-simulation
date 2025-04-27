function [dataflow, prmflow, status] = reconnode_HelicalBackprojection(dataflow, prmflow, status)
% recon node, Helical BP (afrer reconnode_BPprepare)
% [dataflow, prmflow, status] = reconnode_HelicalBackprojection(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_helicalbackprojectionprepare(dataflow, prmflow, status);
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
    carrynode = status.currentjob.carrynode;
end

% main
if pipeline_onoff
    [dataflow.pipepool.(carrynode)(1).data, dataflow.buffer.(nodename)] = ...
        HelicalBPKernelfuntion(dataflow.pipepool.(carrynode)(1).data, dataflow.pipepool.(nodename)(1).data, ...
        dataflow.buffer.(nodename), dataflow.pipepool.(nodename)(1), dataflow.pipepool.(nextnode)(1));
else
    dataflow = HelicalBPKernelfuntion(dataflow, dataflow);
end

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done
1;

% Kernel funtion
    function [dataOut, buffer] = HelicalBPKernelfuntion(dataOut, dataIn, buffer, currpool, nextpool)
        % The anonymous function is static
        debug = [];

        % parameters
        recon = prmflow.recon;
        Npixel = recon.Npixel;
        Npixel_up = recon.Npixel_up;
        Nslice = recon.Nslice;
        isXupsampling = recon.upsampling;
        Xupsampgamma = recon.upsampgamma;
        delta_d_up = recon.delta_d_up/recon.SID;
        midchannel_up = recon.midchannel_up;
        imagesperpitch = recon.imagesperpitch;
        ConeWeightScale = recon.ConeWeightScale;
        Nviewprot = recon.Nviewprot;
%         Zrange = [0 recon.Zviewend];  % I know Zviewstart=0
        ViewRange = [0 recon.Nviewact-1];
        Zscale = recon.imageincrement/recon.delta_z;
        imageZgrid = recon.imageZgrid;
        Zviewshift = recon.Zviewshift;
        % Z-upsampling
        Zupsamp = recon.Zupsamp;
        % recon method
        ispiline = strcmpi(prmflow.recon.method, 'helicalpiline');
        % is concyclic
        concyclic = recon.concyclic;
        % GPU
        GPU_onoff = ~isempty(status.GPUinfo);
%         sigma_z = (Nslice-1)/2/Cd/Zscale;
%         ConeWeightScale_Cz = ConeWeightScale*Cd*Zscale;
        
        reconcenter = recon.center;
        XY = recon.XY;
        Sxy = recon.activeXY;     % Sxy
        NactiveXY = recon.NactiveXY;  % Nxy
        imagesize = recon.imagesize;
        imagesize2 = imagesize(2)*imagesize(1);
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
            % BV
            BV = iteration.BV;
            % noise enhance
            NoiseEnhance = iteration.NoiseEnhance;
        end

        % data classes
        if GPU_onoff
            dataclass_single = ones(1, 'single', 'gpuArray');
            dataclass_raw = gpuArray(ones(1, 'like', dataIn.rawdata));
        else
            dataclass_single = ones(1, 'single');
            dataclass_raw = ones(1, 'like', dataIn.rawdata);
        end
        if iteration_onoff
            dataclass_img = dataclass_raw*(1+1i);
        else
            dataclass_img = dataclass_raw;
        end

        % filters
        basicfilter = recon.filter.basicfilter;
        fftlength = length(basicfilter);

        % ini image
        if ~isfield(dataOut, 'image')
            dataOut.image = zeros(imagesize(1)*imagesize(2), 0, 'like', dataclass_img);
        end
        % ini imagehead
        if ~isfield(dataOut, 'imagehead') || isempty(fieldnames(dataOut.imagehead))
            if pipeline_onoff
                dataOut.imagehead = initialimagehead(nextpool.poolsize);
            else
                dataOut.imagehead = initialimagehead(recon.Nimage);
            end
        end
        if iteration_onoff && ~isfield(dataOut.imagehead, 'Residual')
            if pipeline_onoff
                dataOut.imagehead.Residual = zeros(iteration.Niterloop, nextpool.poolsize, 'single');
            else
                dataOut.imagehead.Residual = zeros(iteration.Niterloop, recon.Nimage, 'single');
            end
        end

        % the console parameters for pipeline on or off
        if pipeline_onoff
            plconsol = status.currentjob.pipeline;
            index_in = poolindex(currpool, plconsol.Index_in);
            index_out = poolindex(nextpool, plconsol.Index_out);
            % Zview_in = plconsol.Do2i;
            % imageindex = index_out - index_out(1) + 1 + nextpool.WritePoint - nextpool.WriteStart;
        else
            index_in = recon.Nviewskip+1 : recon.Nview;
            index_out = 1 : recon.Nimage;
            % Zview_in = 0;
            % imageindex = 1 : recon.Nimage;
        end
        NviewIn = length(index_in);
        Nimageout = length(index_out);
        Zview_in = recon.viewread;
        imageindex = recon.imagewritten + (1 : Nimageout);
        Zgridout = imageZgrid(imageindex);

        % record view read and image written
        if pipeline_onoff
            prmflow.recon.viewread = prmflow.recon.viewread + NviewIn;
            prmflow.recon.imagewritten = prmflow.recon.imagewritten + plconsol.writenumber;
        else
%             prmflow.recon.viewread = prmflow.recon.viewread + NviewIn;
%             prmflow.recon.imagewritten = prmflow.recon.imagewritten + Nimageout;
        end

        % couch direction
        if recon.couchdirection<0
            % backward couch
            index_slice = 1 : Nslice*Npixel;
        else
            % forward couch
            index_slice = reshape(fliplr(reshape(1 : Nslice*Npixel, Npixel, [])), 1, []);
        end
        
        % ini data
        data_f = zeros(fftlength, Nslice*NviewIn, 'like', dataclass_raw);

        % get FBP data
        if GPU_onoff
            data0 = gpuArray(dataIn.rawdata(index_slice, index_in));
        else
            data0 = dataIn.rawdata(index_slice, index_in);
        end
        AngleEncoder0 = dataIn.rawhead.Angle_encoder(index_in);

        % X-upsampling
        if isXupsampling
            data_f(1:Npixel_up, :) = doubleup(reshape(data0, Npixel, []), Xupsampgamma);
        else
            data_f(1:Npixel_up, :) = reshape(data0, Npixel_up, []);
        end
        % view angle
        viewangle = single(double(AngleEncoder0).*anglepercode + viewangleshift + pi/2);
        % sin cos
        costheta = cos(viewangle);
        sintheta = sin(viewangle);

        % branch of right part
        if iteration_onoff
            eta_C = reconcenter(1).*sintheta - reconcenter(2).*costheta;
            channelstart = [floor(recon.midchannel + (-recon.effFOV/2 + eta_C)./recon.delta_d); ...
                          floor(recon.midchannel - (-recon.effFOV/2 + eta_C)./recon.delta_d)];
            if pipeline_onoff
                spaceshift = dataflow.buffer.(nodename).spaceshift;
                for isub = 1:iteration.viewsubspace
                    % the branch node
                    branchnode = plconsol.branchoutput(isub).nextnode;
                    branchpoolindex = plconsol.branchoutput(isub).poolindex;
                    branchcarrynode = plconsol.branchoutput(isub).carrynode;
                    branchcarryindex = plconsol.branchoutput(isub).carryindex;
                    % branch subviewshift
                    subviewshift = iteration.subviewshift(isub);
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
                    dataflow.pipepool.(branchcarrynode)(branchcarryindex).data.rawdata(:, index_branch) = ...
                        data0(:, index_subview);
                    dataflow.pipepool.(branchcarrynode)(branchcarryindex).data.rawhead.Angle_encoder(index_branch) = ...
                        AngleEncoder0(index_subview);
                    % channelstart
                    dataflow.pipepool.(branchcarrynode)(branchcarryindex).data.rawhead.channelstart(:, index_branch) = ...
                        channelstart(:, index_subview);
                end
                % update the spaceshift
                dataflow.buffer.(nodename).spaceshift = mod(spaceshift + (ceil(NviewIn/iteration.viewsubspace)* ...
                    iteration.viewsubspace - NviewIn), iteration.viewsubspace);
            else
                for ifield = fieldnames(dataOut.rawhead)'
                    dataOut.rawhead.(ifield{1}) = dataOut.rawhead.(ifield{1})(:, index_in);
                end
                dataOut.rawhead.channelstart = channelstart;
                dataOut.rawdata = dataOut.rawdata(:, index_in);
            end
        end

        % BP prepare
        % ZviewIndex
        ZviewIdx = cast(((0:NviewIn-1)+Zview_in), 'like', dataclass_single);
        Nviewprot_GPU = cast(Nviewprot, 'like', dataclass_single);
        ViewRange = cast(ViewRange, 'like', dataclass_single);

        % Filter, inner node function
        if ~recon.filtered
            % isreal
            isreal_flag = isreal(data_f);
            if isreal_flag
                data_f = ifft(fft(data_f).*basicfilter, 'symmetric');
            else
                data_f = ifft(fft(data_f).*basicfilter);
            end
        end
        data_f = reshape(data_f, fftlength, Nslice, NviewIn);
        
        % copy some variable to GPU
        Nslice_GPU = gpuArray(single(Nslice));
        Cd = cast(imagesperpitch, 'like', dataclass_single);
        costheta = cast(costheta, 'like', dataclass_single);
        sintheta = cast(sintheta, 'like', dataclass_single);
        % BP normalize
        BPnormalize = cast(pi/Nviewprot, 'like', dataclass_single);

        % ini image data to BP
        image_fix = zeros(NactiveXY, Nimageout, 'like', dataclass_raw);
        % loop the views
        for iview = 1:NviewIn
        %for iview = 1:125
        % for iview = 1:20:NviewIn
            % X-Y to Zeta-Eta
            Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
            Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);

            % interp target on Z
            Deta = sqrt(1-Eta.^2);
            Phi = asin(Eta)./(pi*2);
            Zv = ZviewIdx(iview);

            if ispiline
                % Pi-line BP
                [Tz, Tchn, Weight] = pilineinterp(Zgridout-Zviewshift, Zv, Eta, Zeta, Deta, Phi, Cd, Nviewprot_GPU, ...
                    Nslice_GPU, Zscale, Zupsamp.Zupsampling, delta_d_up, midchannel_up, Nimageout, concyclic);
            else
                % normal BP (not support concyclic yet)
                [Tz, Tchn, Weight] = helicalbpinterp(Zgridout-Zviewshift, Zv, Eta, Zeta, Deta, Phi, Cd, Nviewprot_GPU, ...
                    Nslice_GPU, Zscale, Zupsamp.Zupsampling, delta_d_up, midchannel_up, Nimageout, ConeWeightScale, ViewRange);
            end
            
            % Z upsampling of BP data
            D = data_f(1:Npixel_up, :, iview)*Zupsamp.ZupMatrix;
            % 2D interp on channels and slices
            data_2 = interp2(D, Tz, Tchn, 'linear', 0);

            % add to image
            image_fix = image_fix + data_2.*Weight;
        end
        % normalize by pi/Nviewprot
        image_fix = image_fix.*BPnormalize;
        % add to output
        if size(dataOut.image, 2) < max(index_out)
            dataOut.image(:, max(index_out)) = 0;
        end
        dataOut.image(Sxy(:), index_out) = dataOut.image(Sxy(:), index_out) + image_fix;

        % image head and TV
        if pipeline_onoff && plconsol.newAvail > 0
            % avialable images index
            index_avail = poolindex(nextpool, plconsol.Index_avail);
            % imageindex_avail = index_avail - nextpool.WriteStart + 1;
            imageindex_avail = recon.availwritten + (1 : plconsol.newAvail);
            % record avialable image written
            prmflow.recon.availwritten = prmflow.recon.availwritten + plconsol.newAvail;
            % head
            dataOut.imagehead = reconimagehead(dataOut.imagehead, recon, imageindex_avail, index_avail);
            % TV
            if iteration_onoff
                % get images to do TV
                index_TV = index_avail;
                Ntv = plconsol.newAvail;
                u = cast(reshape(real(dataOut.image(:, index_TV)), imagesize(2), imagesize(1), Ntv), 'like', dataclass_img);
                u = u.*(1 + 1i);
                % TV
                if index_avail(1) > 1                  
                    u = cat(3, buffer.ub, u);
                    [u, b] = BregmanTV3D_helicaliter(u, BV.mu_BV0, BV.Cl, [], BV.Crange, BV.maxNiter, ...
                        BV.tol2, BV.Zbalance);
                    mu_BV = BV.mu_BV0./(abs(real(u) - imag(u)).*BV.mu_adap + 1);
                    u = BregmanTV3D_helicaliter(u, mu_BV, BV.Cl, b, BV.Crange, BV.maxNiter, BV.tol, BV.Zbalance);
                else
                    [u, b] = BregmanTV3D_axialiter(u, BV.mu_BV0, BV.Cl, [], BV.Crange, BV.maxNiter, BV.tol2, ...
                        BV.Zbalance);
                    mu_BV = BV.mu_BV0./(abs(real(u) - imag(u)).*BV.mu_adap + 1);
                    u = BregmanTV3D_axialiter(u, mu_BV, BV.Cl, b, BV.Crange, BV.maxNiter, BV.tol, BV.Zbalance);
                end
                % back up boundary of u
                buffer.ub = u(:,:,end-1);
                % noise enhance
                u = enhancemix_iter(u, NoiseEnhance.wlevel, NoiseEnhance.mixwidth, NoiseEnhance.mixalpha, NoiseEnhance.mixbeta);
                % or, u = real(u) + (real(u)-imag(u)).*(-1/2+1/2i);
                % output
                u = reshape(u, imagesize2, []);
                if index_avail(1) > 1
                    if plconsol.isshotend
                        dataOut.image(:, index_avail) = u(:, 2:end);
                    else
                        dataOut.image(:, index_avail(1:end-1)) = u(:, 2:end-1);
                    end
                else
                    if plconsol.isshotend
                        dataOut.image(:, index_avail) = u;
                    else
                        dataOut.image(:, index_avail(1:end-1)) = u(:, 1:end-1);
                    end
                end
                % fix newAvail
                if ~plconsol.isshotend
                    status.currentjob.pipeline.newAvail = status.currentjob.pipeline.newAvail-1;
                end
            end
            
        elseif ~pipeline_onoff
            % head
            dataOut.imagehead = reconimagehead(dataOut.imagehead, recon, imageindex);
            % TV for iteration
            if iteration_onoff
                u = cast(reshape(dataOut.image, imagesize(2), imagesize(1), Nimageout), 'like', dataclass_img);
                [Gu, b] = BregmanTV3D_axialiter(u, BV.mu_BV0, BV.Cl, [], [], BV.Crange, BV.maxNiter, BV.tol2, ...
                        BV.Zbalance);
                mu_BV = BV.mu_BV0./(abs(u - Gu).*BV.mu_adap + 1);
                Gu = BregmanTV3D_axialiter(real(u), mu_BV, BV.Cl, Gu, b, BV.Crange, BV.maxNiter, BV.tol, BV.Zbalance);
                u = (real(u) + Gu)./2 + 1.0i.*(real(u) - Gu)./2;
                dataOut.image = gather(reshape(u, [], Nimageout));
            end
        end
        

        % jobdone
        if ~pipeline_onoff
            status.jobdone = true;
        end
        % BP Kernelfuntion END
    end

end


function [Tz, Tchn, Weight] = pilineinterp(Zgrid_f, Zv, Eta, Zeta, Deta, Phi, Cd, Nprot, Nslice, Zscale, Zupsampling, ...
    delta_d, midchannel, Nimage, concyclic)
% pi-line BP

if nargin < 15
    concyclic = false;
    % concyclic-mode 
end

if ~concyclic
    % Tz
    Tz = (Zgrid_f - (Zv/Nprot - Phi).*Cd)./(Deta+Zeta);
    % in Pi?
    PiC = (Tz.*Deta./Cd - Phi).*4;
else
    Tz = (Zgrid_f - (Zv/Nprot - Phi).*Cd).*Deta./(Deta+Zeta);
    PiC = (Tz./Cd - Phi).*4;
end
Weight = PiC>=-1 & PiC<1;

% Z scale
Tz = Tz.*Zscale + (Nslice+1)/2;
% extrap (for big pitch)
% Tz(Tz<1) = 1;  Tz(Tz>Nslice) = Nslice;
Tz = Tz.*(Tz>=1 & Tz<=Nslice) + (Tz<1).*1.0 + (Tz>Nslice).*Nslice;
% shift Tz to the Z-upsampled position
Tz = (Tz - 1).*Zupsampling + 1;

% interp target on Eta
Tchn = repmat(Eta./delta_d + midchannel, 1, Nimage);

end




function [dataflow, prmflow, status] = BPpriostep(dataflow, prmflow, status)
% Helical BP priostep

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
nextnode = status.pipeline.(nodename).nextnode;
pipeprm = nodeprm.pipeline;

imagenumber = double(prmflow.recon.imagenumber);
% Nview = double(prmflow.recon.Nview);
Nviewskip = double(prmflow.recon.Nviewskip);
Nviewact = double(prmflow.recon.Nviewact);
switch lower(prmflow.recon.method)
    case {'helical', 'helical3d'}
        % normal Helical
        viewbyimages = prmflow.recon.viewbyimages_full;
    case 'helicalpiline'
        % Helical pi-line
        viewbyimages = prmflow.recon.viewbyimages_pi;
    otherwise
        % error
        error('Unknown reconstruction method %s for helical!', prmflow.recon.method);
end

% iteration
if isfield(prmflow.recon, 'iteration_onoff')
    iteration_onoff = prmflow.recon.iteration_onoff;
else
    iteration_onoff = false;
end

% ini the flag to run the poststep
status.currentjob.torunpoststep = true;

% prio step #1, the prio steps of the pipepools
% pass due to nextpool is stucked
if dataflow.pipepool.(nextnode)(1).WriteStuck
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 6;
    % done and keep waking
    return;
end

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

% I know the nodetype is H-H.1.G

% check shot start 1
status.currentjob.pipeline.isshotstart = dataflow.pipepool.(nextnode)(1).isshotstart;
% close it
dataflow.pipepool.(nextnode)(1).isshotstart = false;
% Isshotstart1 = dataflow.pipepool.(nodename)(1).ReadPoint == dataflow.pipepool.(nodename)(1).ReadStart;

if status.currentjob.pipeline.isshotstart
    % ini next pool
    dataflow.pipepool.(nextnode)(1).ReadStart = 1;
    dataflow.pipepool.(nextnode)(1).ReadEnd = imagenumber;
    dataflow.pipepool.(nextnode)(1).WriteStart = 1;
    dataflow.pipepool.(nextnode)(1).WriteEnd = imagenumber;
    dataflow.pipepool.(nextnode)(1).ReadPoint = 1;
    dataflow.pipepool.(nextnode)(1).WritePoint = 1;
    % close the isshotstart
    dataflow.pipepool.(nextnode)(1).isshotstart = false;
    % view skip
    dataflow.pipepool.(nodename)(1).ReadPoint = dataflow.pipepool.(nodename)(1).ReadStart + Nviewskip;
    dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd - Nviewskip;
    % Note: we moved the ReadPoint and ReadEnd but not the ReadStart
    % branches
    if iteration_onoff
        viewsubspace = prmflow.iteration.viewsubspace;
        subviewshift = prmflow.iteration.subviewshift;
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
                % skip to close the isshotstart
                dataflow.pipepool.(nextnode)(1).isshotstart = true;
                % to pass the kernel function
                status.currentjob.topass = true;
                % withdraw the initial of next pool
                if status.currentjob.pipeline.isshotstart
                    dataflow.pipepool.(nextnode)(1).isshotstart = true;
                    dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
                end
                return;
            end
            Nview_ii = double(floor((prmflow.recon.Nviewact - subviewshift(ii)) / viewsubspace) + 1);
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadEnd = Nview_ii;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteStart = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WriteEnd = Nview_ii;
            dataflow.pipepool.(branchnode)(branchpoolindex).ReadPoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).WritePoint = 1;
            dataflow.pipepool.(branchnode)(branchpoolindex).AvailPoint = 0;
        end
        1;
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

% check shot start 2
% status.currentjob.pipeline.isshotstart = dataflow.pipepool.(nodename)(1).ReadPoint == ...
%     dataflow.pipepool.(nodename)(1).ReadStart + Nviewskip;

% check shot end 1 (the input data has reached the end)
Isshotend1 = dataflow.pipepool.(nodename)(1).WritePoint == dataflow.pipepool.(nodename)(1).WriteEnd + 1;

% minlimit/maxlimit
status.currentjob.pipeline.minlimit = pipeprm.inputminlimit;
status.currentjob.pipeline.maxlimit = pipeprm.inputmaxlimit;

% avail inputs
currAvailNumber = min(dataflow.pipepool.(nodename)(1).AvailPoint, dataflow.pipepool.(nodename)(1).ReadEnd) ...
    - dataflow.pipepool.(nodename)(1).ReadPoint + 1;
% I know, we have moved the ReadEnd therefore the AvailPoint could run over the ReadEnd.
n = min(pipeprm.inputmaxlimit, currAvailNumber);

% space left in nextpool
nextleft = dataflow.pipepool.(nextnode)(1).poolsize - dataflow.pipepool.(nextnode)(1).WritePoint + 1;

% m (active images)
AvailV0 = [dataflow.pipepool.(nodename)(1).ReadPoint dataflow.pipepool.(nodename)(1).AvailPoint] - ...
    dataflow.pipepool.(nodename)(1).ReadStart - Nviewskip + 1;
AvailV = [0  dataflow.pipepool.(nodename)(1).AvailPoint - dataflow.pipepool.(nodename)(1).ReadPoint] + ...
    prmflow.recon.viewread + 1;
% I know the reading was starting from dataflow.pipepool.(nodename)(1).ReadStart+Nviewskip
s = (viewbyimages(2, :) >= AvailV(1)) & (viewbyimages(1, :) <= AvailV(2));
m = sum(s);
if m > nextleft
    m = nextleft;
    n = viewbyimages(1, find(s, 1) + nextleft) - AvailV(1);
    if n < status.currentjob.pipeline.minlimit || nextleft <= 0
        % stucking
        status.currentjob.pipeline.readnumber = 0;
        status.currentjob.pipeline.writenumber = 0;
        status.currentjob.pipeline.newAvail = 0;
        status.jobdone = 6;
        % to pass the kernel function
        status.currentjob.topass = true;
        % withdraw the initial of next pool
        if status.currentjob.pipeline.isshotstart
            dataflow.pipepool.(nextnode)(1).isshotstart = true;
            dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
        end
        return;
    end
elseif n < status.currentjob.pipeline.minlimit && ~Isshotend1
    % not enough input views
    status.currentjob.pipeline.readnumber = 0;
    status.currentjob.pipeline.writenumber = 0;
    status.currentjob.pipeline.newAvail = 0;
    status.jobdone = 3;
    % to pass the kernel function
    status.currentjob.topass = true;
    % withdraw the initial of next pool
    if status.currentjob.pipeline.isshotstart
        dataflow.pipepool.(nextnode)(1).isshotstart = true;
        dataflow.pipepool.(nodename)(1).ReadEnd = dataflow.pipepool.(nodename)(1).ReadEnd + Nviewskip;
    end
    return;
end
if n < currAvailNumber
    % I know nextleft>0
    if currAvailNumber - n >= pipeprm.inputminlimit || Isshotend1
        % partly done
        status.jobdone = 2;
    else
        status.jobdone = 1;
    end
else
    % done
    status.jobdone = 1;
end

% check shot end 2 (will output the shot end)
status.currentjob.pipeline.isshotend = Isshotend1 && n == currAvailNumber;

% Index_in, Index_out
status.currentjob.pipeline.Index_in = ...
    [dataflow.pipepool.(nodename)(1).ReadPoint dataflow.pipepool.(nodename)(1).ReadPoint+n-1];
status.currentjob.pipeline.Index_out = ...
    [dataflow.pipepool.(nextnode)(1).WritePoint dataflow.pipepool.(nextnode)(1).WritePoint+m-1];
% Do2i
% status.currentjob.pipeline.Do2i = ...
%     dataflow.pipepool.(nodename)(1).ReadPoint - dataflow.pipepool.(nodename)(1).ReadStart - Nviewskip;
% % I know the Do2i will be used to calculate the Zview, not what it was.
status.currentjob.pipeline.Do2i = 0;
% Nexpand (has been count)
status.currentjob.pipeline.Nexpand = 0;

% readnumber
status.currentjob.pipeline.readnumber = n;

% writenumber
s = (viewbyimages(2, :) >= AvailV(1)) & (min(viewbyimages(2, :), Nviewact) <= AvailV(1)+n-1);
writenumber = sum(s);
status.currentjob.pipeline.writenumber = writenumber;
% newAvail
% status.currentjob.pipeline.newAvail = writenumber;
Navail = dataflow.pipepool.(nextnode)(1).WritePoint - dataflow.pipepool.(nextnode)(1).AvailPoint - 1 + writenumber;
if Navail < pipeprm.outputminlimit && ~status.currentjob.pipeline.isshotend
    status.currentjob.pipeline.newAvail = 0;
else
    status.currentjob.pipeline.newAvail = Navail;
end
% index_avail
status.currentjob.pipeline.Index_avail = ...
    dataflow.pipepool.(nextnode)(1).AvailPoint + [1 status.currentjob.pipeline.newAvail];
    
% fix jobdone
if status.currentjob.pipeline.newAvail == 0
    if status.jobdone == 1
        status.jobdone = 4;
    elseif status.jobdone == 2
        status.jobdone = 7;
    end
end

% branchs in status.currentjob
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


