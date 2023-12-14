function [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status)
% recon node, Axial 3D BP (afrer reconnode_BPprepare)
% [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status);

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

% BP parameters
nodeprm = prmflow.pipe.(status.nodename);
% Filter?
isFilter = isfield(nodeprm, 'Filter');
% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end
% echo onoff
if isfield(status, 'echo_onoff')
    echo_onoff = status.echo_onoff;
else
    echo_onoff = false;
end

% recon parameters
recon = prmflow.recon;
% SID = recon.SID;
Nshot = recon.Nshot;
% FOV = recon.FOV;
Nviewprot = recon.Nviewprot;
% startviewangle = recon.startviewangle;
imagesize = recon.imagesize;
% midchannel = recon.midchannel;
% delta_d = recon.delta_d/SID;
Npixel = recon.Npixel;
Nslice = recon.Nslice;
Neighb = recon.Neighb;
Nextslice = recon.Nextslice;
Nimage = recon.Nimage;
% reconcenter = recon.center;
couchdirection = recon.couchdirection;
viewblock = recon.viewblock;
Nviewblock = ceil(Nviewprot/viewblock);
XY = recon.XY;
Sxy = recon.activeXY;
NactiveXY = recon.NactiveXY;

% is upsampling?
if recon.upsampling || recon.upsampled
    Npixel_up = recon.Npixel_up;
    midchannel = recon.midchannel_up;
    delta_d = recon.delta_d_up/recon.SID;
else
    Npixel_up = Npixel;
    delta_d = recon.delta_d/recon.SID;
    midchannel = recon.midchannel;
end
% warning when the filter node is suspected in wrong order
if recon.upsampling && ~isFilter
    warning('The upsampling shall be done before the filter!');
end
upsampgamma = recon.upsampgamma;

% is filtering
if isFilter
    filter = recon.filter;
    if isstruct(filter)
        filter = filter.basicfilter;
    end
    filtlen = length(filter);
else
    filter = [];
    filtlen = [];
end

% viewangle
viewangle = dataflow.rawhead.viewangle;
% rawdata reshape
dataflow.rawdata = reshape(dataflow.rawdata, [], Nslice, Nviewprot, Nshot);
if couchdirection>0
    % forward couch
    dataflow.rawdata = flip(dataflow.rawdata, 2);
end

% ini image
dataflow.image = zeros(imagesize(1)*imagesize(2), Nimage, 'single');

% channel interp  
Chninterp.delta_d = delta_d;
Chninterp.midchannel = midchannel;
Chninterp.channelindex = single(1:Npixel_up)';

% Zinterp (it was prepared in reconnode_BPprepare) 
Zinterp = recon.Zinterp;
Nfillslice = Nslice + Zinterp.Nfill0*2;

% Z-upsampling
Zupsamp = recon.Zupsamp;
Nslice_up = (Nfillslice-1)*Zupsamp.Zupsampling+1;
if Zupsamp.Zupsampling==0
    error('The Zupsampling shall not be 0!');
end

% to GPU
if GPUonoff
    [XY, Chninterp, Zinterp, upsampgamma, Zupsamp, filter, filtlen, Npixel, Npixel_up] = putinGPU(...
     XY, Chninterp, Zinterp, upsampgamma, Zupsamp, filter, filtlen, Npixel, Npixel_up);
end

for ishot = 1:Nshot
%     tic;
    % image index and shotflag
    [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot);
    % three shotflags
    switch shotflag
        case 1
            % first shot
            Nactslice = Nslice + Neighb;
            Ztarget = -Nslice/2+1 : Nslice/2+Neighb;
            ZupMatrix = Zupsamp.ZupMatrix_1;
        case 2
            % last shot
            Nactslice = Nslice + Neighb;
            Ztarget = -Nslice/2-Neighb+1 : Nslice/2;
            ZupMatrix = Zupsamp.ZupMatrix_end;
        case 3
            % only one shot
            Nactslice = Nslice;
            Ztarget = -Nslice/2+1 : Nslice/2;
            ZupMatrix = Zupsamp.ZupMatrix_1 + Zupsamp.ZupMatrix_end - Zupsamp.ZupMatrix_mid;
        otherwise
            % middle shots
            Nactslice = Nslice + Neighb*2;
            Ztarget = -Nactslice/2+1 : Nactslice/2;
            ZupMatrix = Zupsamp.ZupMatrix_mid;
    end
    % ini
    img_shot = zeros(NactiveXY, Nactslice, 'single');
    Ztarget = repmat(single(Ztarget), NactiveXY, 1);
    if GPUonoff
        [img_shot, Ztarget] = putinGPU(img_shot, Ztarget);
    end
    Nfill0 = Zinterp.Nfill0;
    % loop the view blocks
    for iblk = 1:Nviewblock
        if iblk<Nviewblock
            viewindex = (1:viewblock) + (iblk-1)*viewblock;
        else  % iblk == Nviewblock
            viewindex = (iblk-1)*viewblock+1 : Nviewprot;
        end
        Nviewperblk = length(viewindex);
        % get data per block
        if GPUonoff
            datablk = gpuArray(dataflow.rawdata(:, :, viewindex, ishot));
            viewangleblk = gpuArray(viewangle(viewindex));
        else
            datablk = dataflow.rawdata(:, :, viewindex, ishot);
            viewangleblk = viewangle(viewindex);
        end
        % up sampling
        if recon.upsampling
            datablk = doubleup(reshape(datablk, Npixel, []), upsampgamma);
        end
        % filter
        if isFilter
            datablk_f = fft(reshape(datablk, Npixel_up, []), filtlen);
            datablk_f = ifft(datablk_f.*filter, 'symmetric');
            datablk = reshape(datablk_f(1:Npixel_up, :), Npixel_up, Nslice, []);
        end
        % Z upsampling
        datablk = reshape(permute(datablk, [1 3 2]), Npixel_up*Nviewperblk, Nslice);
        datablk = datablk*ZupMatrix;
        datablk = permute(reshape(datablk, Npixel_up, Nviewperblk, Nslice_up), [1 3 2]);

        costheta = cos(viewangleblk);
        sintheta = sin(viewangleblk);
        % loop the views
        for iview = 1:Nviewperblk
            % X-Y to Zeta-Eta
            Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
            Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);

            % interp target on Eta
            Tchn = repmat(Eta./Chninterp.delta_d + Chninterp.midchannel, 1, Nactslice);
            
            % interp target on Z
            Tz = interp3(Zinterp.Zeta, Zinterp.Eta, Zinterp.zz, Zinterp.t, ...
                repmat(Zeta, 1, Nactslice), repmat(Eta, 1, Nactslice), Ztarget);
            switch shotflag
                case 1
                    % first shot
                    % Tz(Tz<Nfill0+1) = Nfill0+1;
                    Tz = Tz.*(Tz>=Nfill0+1) + (Tz<Nfill0+1).*(Nfill0+1);
                case 2
                    % last shot
                    % Tz(Tz>Nfill0+Nslice) = Nfill0+Nslice;
                    Tz = Tz.*(Tz<=Nfill0+Nslice) + (Tz>Nfill0+Nslice).*(Nfill0+Nslice);
                case 3
                    % only one shot
                    % Tz(Tz<Nfill0) = Nfill0; Tz(Tz>Nfill0+Nslice) = Nslice;
                    Tz = Tz.*((Tz>=Nfill0+1) & (Tz<=Nfill0+Nslice)) + (Tz<Nfill0+1).*(Nfill0+1) + (Tz>Nfill0+Nslice).*(Nfill0+Nslice);
                otherwise
                    % do nothing
                    1;
            end
            % I know Nfill0 = 2.

            Tz = (Tz-1).*Zupsamp.Zupsampling + 1;
            img_shot = img_shot + interp2(datablk(:,:, iview), Tz, Tchn, 'linear', 0);
   
        end

        % echo '.'
        if echo_onoff, fprintf('.'); end
    end
    % get img
    dataflow.image(Sxy ,imageindex) = dataflow.image(Sxy ,imageindex) + gather(img_shot);
    
end
% normalize by Nviewprot, (pi/2 due to the filter)
dataflow.image = reshape(dataflow.image, imagesize(2), imagesize(1), Nimage).*(pi/Nviewprot/2);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot)

imageindex = (1:Nextslice) + ((ishot-1)*Nslice-(Nextslice-Nslice)/2);
if Nshot==1
    imageindex(1:(Nextslice-Nslice)/2) = [];
    imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
    shotflag = 3;
elseif ishot==1
    imageindex(1:(Nextslice-Nslice)/2) = [];
    shotflag = 1;
elseif ishot==Nshot
    imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
    shotflag = 2;
else
    shotflag = 0;
end

end