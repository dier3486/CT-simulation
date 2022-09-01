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
SID = prmflow.recon.SID;
Nshot = prmflow.recon.Nshot;
% FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
startviewangle = prmflow.recon.startviewangle;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d/SID;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Neighb = prmflow.recon.Neighb;
Nextslice = prmflow.recon.Nextslice;
Nimage = prmflow.recon.Nimage;
% reconcenter = prmflow.recon.center;
couchdirection = prmflow.recon.couchdirection;
viewblock = prmflow.recon.viewblock;
Nviewblock = ceil(Nviewprot/viewblock);
XY = prmflow.recon.XY;
Sxy = prmflow.recon.activeXY;
NactiveXY = prmflow.recon.NactiveXY;
viewangle = prmflow.recon.viewangle;
upsampling = prmflow.recon.upsampling;
upsampgamma = prmflow.recon.upsampgamma;
if isFilter
    filter = prmflow.recon.filter;
else
    filter = [];
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

% ini image
dataflow.image = zeros(imagesize*imagesize, Nimage, 'single');

% channel interp  
Chninterp.delta_d = delta_d;
Chninterp.midchannel = midchannel;
Chninterp.channelindex = single(1:Npixel_up)';

% Zinterp (it was prepared in reconnode_BPprepare) 
Zinterp = prmflow.recon.Zinterp;
Nfillslice = Nslice + Zinterp.Nfill0*2;

% to GPU
if GPUonoff
    [XY, Chninterp, Zinterp, upsampgamma, filter, Npixel, Npixel_up] = ...
        putinGPU(XY, Chninterp, Zinterp, upsampgamma, filter, Npixel, Npixel_up);
end
filtlen = length(filter);

for ishot = 1:Nshot
%     tic;
    % image index and shotflag
    [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot, couchdirection);
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
    % loop the view blocks
    for iblk = 1:Nviewblock
        if iblk<Nviewblock
            viewindex = (1:viewblock) + (iblk-1)*viewblock;
        else  % iblk == Nviewblock
            viewindex = (iblk-1)*viewblock+1 : Nviewprot;
        end       
        % get data per block
        if GPUonoff
            datablk = gpuArray(dataflow.rawdata(:, :, viewindex, ishot));
            viewangleblk = gpuArray(viewangle(viewindex) + startviewangle(ishot));
        else
            datablk = dataflow.rawdata(:, :, viewindex, ishot);
            viewangleblk = viewangle(viewindex) + startviewangle(ishot);
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
        img_shot = backproject3D(img_shot, datablk, viewangleblk, XY, Chninterp, Zinterp, shotflag, interpmethod, ...
            index_img, Ztarget, databuff1, databuff2);
        % echo '.'
        if echo_onoff, fprintf('.'); end
    end
    % get img
    dataflow.image(Sxy ,imageindex) = dataflow.image(Sxy ,imageindex) + gather(img_shot);
    
end
% normalize by Nviewprot, (pi/2 due to the filter)
dataflow.image = reshape(dataflow.image, imagesize, imagesize, Nimage).*(pi/Nviewprot/2);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot, couchdirection)

if couchdirection<0
    % backward couch
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
else
    % forward couch
    imageindex = (1:Nextslice) + ((Nshot-ishot)*Nslice-(Nextslice-Nslice)/2);
    if Nshot==1
        imageindex(1:(Nextslice-Nslice)/2) = [];
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        shotflag = 3;
    elseif ishot==Nshot
        imageindex(1:(Nextslice-Nslice)/2) = [];
        shotflag = 1;
    elseif ishot==1
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        shotflag = 2;
    else
        shotflag = 0;
    end
end

end