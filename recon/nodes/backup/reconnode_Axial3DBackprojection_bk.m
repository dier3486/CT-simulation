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
% interp method (linear or 4points)
if isfield(nodeprm, 'interp')
    interpmethod = nodeprm.interp;
else
    interpmethod = '4points';
end
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
SID = recon.SID;
Nshot = recon.Nshot;
% FOV = recon.FOV;
Nviewprot = recon.Nviewprot;
% startviewangle = recon.startviewangle;
imagesize = recon.imagesize;
midchannel = recon.midchannel;
delta_d = recon.delta_d/SID;
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
viewangle = recon.viewangle;
upsampling = recon.upsampling;
upsampgamma = recon.upsampgamma;
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
if upsampling
    Npixel_up = Npixel*2;
    delta_d = delta_d/2;
    midchannel = midchannel*2-1;
    if ~isFilter
        % It is a mistake!
        warning('The upsampling shall be done before the filter!');
    end
else
    Npixel_up = Npixel;
end

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);
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
Zupsampling = Zupsamp.Zupsampling;
if Zupsampling>1
    ZupMatrix = Zupsamp.ZupMatrix;
else
    ZupMatrix = [];
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
        case 2
            % last shot
            Nactslice = Nslice + Neighb;
            Ztarget = -Nslice/2-Neighb+1 : Nslice/2;
        case 3
            % only one shot
            Nactslice = Nslice;
            Ztarget = -Nslice/2+1 : Nslice/2;
        otherwise
            % middle shots
            Nactslice = Nslice + Neighb*2;
            Ztarget = -Nactslice/2+1 : Nactslice/2;
    end
    % ini
    img_shot = zeros(NactiveXY, Nactslice, 'single');
    index_img = repmat(single(1:NactiveXY)', 1, Nactslice);
    Ztarget = repmat(single(Ztarget), NactiveXY, 1);
    if GPUonoff
        [img_shot, index_img, Ztarget] = putinGPU(img_shot, index_img, Ztarget);
        databuff1 = zeros(NactiveXY, Nfillslice, 'single', 'gpuArray');
        databuff2 = zeros(NactiveXY, Nactslice, 'single', 'gpuArray');
    else
        databuff1 = [];
        databuff2 = [];
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
        if upsampling
            datablk = doubleup(reshape(datablk, Npixel, []), upsampgamma);
        end
        % filter
        if isFilter
            datablk_f = fft(reshape(datablk, Npixel_up, []), filtlen);
            datablk_f = ifft(datablk_f.*filter, 'symmetric');
            datablk = reshape(datablk_f(1:Npixel_up, :), Npixel_up, Nslice, []);
        end
        % 3D BP
%         img_shot = backproject3D(img_shot, datablk, viewangleblk, XY, Chninterp, Zinterp, shotflag, interpmethod, ...
%             index_img, Ztarget, databuff1, databuff2);

        costheta = cos(viewangleblk);
        sintheta = sin(viewangleblk);
        % loop the views
        for iview = 1:Nviewperblk
            % X-Y to Zeta-Eta
            Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
            Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);
            % interp on channel direction
            t_chn = Eta./Chninterp.delta_d + Chninterp.midchannel;
            databuff1(:, Nfill0+1:end-Nfill0) = interp1(Chninterp.channelindex, datablk(:, :, iview), t_chn(:), 'linear', 0);
            % start and last shot fillup
            switch shotflag
                case 1
                    % first shot
                    databuff1(:, 1:Nfill0) = repmat(databuff1(:, Nfill0+1), 1, Nfill0);
                case 2
                    % last shot
                    databuff1(:, end-Nfill0+1:end) = repmat(databuff1(:, end-Nfill0), 1, Nfill0);
                case 3
                    % only one shot
                    databuff1(:, 1:Nfill0) = repmat(databuff1(: ,Nfill0+1), 1, Nfill0);
                    databuff1(:, end-Nfill0+1:end) = repmat(databuff1(:, end-Nfill0), 1, Nfill0);
                otherwise
                    % do nothing
                    1;
            end

            % Zinterp index
            Tz = interp3(Zinterp.Zeta, Zinterp.Eta, Zinterp.zz, Zinterp.t, ...
                repmat(Zeta, 1, Nactslice), repmat(Eta, 1, Nactslice), Ztarget);

            if Zupsamp.Zupsampling > 0
            else
            end
            % interp on Z
            switch lower(interpmethod)
                case 'linear'
                    databuff2 = interp2(databuff1, Tz, index_img, 'linear', 0);
                case '4points'
%                     Tz_floor = floor(Tz(:));
%                     s_odd = mod(Tz_floor, 2);
%                     alpha = Tz(:) - Tz_floor;
%                     % beta & gamma
%                     beta = interp1(Zinterp.fourpointindex, Zinterp.fourpoint, alpha.*Zinterp.Nfourp+1);

                    [t_odd, t_even, Gamma] = omiga4interp(Tz, Zupsamp.Cgamma);
%                     % 4-points-interp
%                     t_odd  = reshape((Tz_floor+s_odd)./2 + beta(:,1).*(1-s_odd) + beta(:,2).*s_odd, imagesize2, []);
%                     t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), imagesize2, []);
                    databuff2 = interp2(databuff1(:, 1:2:end), t_odd, index_img, 'linear', 0)./2 + ...
                        interp2(databuff1(:, 2:2:end), t_even, index_img, 'linear', 0)./2;
                    databuff2 = databuff2 + interp2(conv2(databuff1, [-1 2 -1], 'same'), Tz, index_img, 'linear', 0) ...
                        .*Gamma./4;
                    % I know Zinterpsamp.convL = [-1 2 -1]
                otherwise
                    error('Unknown interplation method %s!', interpmethod);
            end
            % add to image
            img_shot = img_shot + databuff2;
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