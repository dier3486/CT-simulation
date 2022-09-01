function [dataflow, prmflow, status] = reconnode_Axial2DBackprojection(dataflow, prmflow, status)
% recon node, Axial 2D BP (afrer reconnode_BPprepare)
% [dataflow, prmflow, status] = reconnode_Axial2DBackprojection(dataflow, prmflow, status);

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
channelindex = single(1:Npixel_up)';

% to GPU
if GPUonoff
    [XY, channelindex, delta_d, midchannel, upsampgamma, filter, Npixel, Npixel_up] = ...
        putinGPU(XY, channelindex, delta_d, midchannel, upsampgamma, filter, Npixel, Npixel_up);
end
filtlen = length(filter);

for ishot = 1 : Nshot
    imageindex = getimageindex(Nslice, Nshot, ishot, couchdirection);
    % ini
    if GPUonoff
        img_shot = zeros(NactiveXY, Nslice, 'single', 'gpuArray');
        databuff = zeros(NactiveXY, Nslice, 'single', 'gpuArray');
    else
        img_shot = zeros(NactiveXY, Nslice, 'single');
        databuff = [];
    end
    % loop the view blocks
    for iblk = 1:Nviewblock
        % viewindex
        if iblk<Nviewblock
            viewindex = (1:viewblock) + (iblk-1)*viewblock;
        else  % iblk == Nviewblock
            viewindex = (iblk-1)*viewblock+1 : Nviewprot;
        end
        Nviewperblk = length(viewindex);
        viewangleblk = viewangle(viewindex) + startviewangle(ishot);
        % get data per block
        if GPUonoff
            datablk = gpuArray(dataflow.rawdata(:, :, viewindex, ishot));
            [viewangleblk, Nviewperblk]  = putinGPU(viewangleblk, Nviewperblk);
        else
            datablk = dataflow.rawdata(:, :, viewindex, ishot);
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
        % 2D BP
        % cos and sin of viewangle
        costheta = cos(viewangleblk);
        sintheta = sin(viewangleblk);
        % loop the views
        for iview = 1:Nviewperblk
            % X-Y to Zeta-Eta
            Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
            % interp on channel direction
            t_chn = Eta./delta_d + midchannel;
            databuff = interp1(channelindex, datablk(:, :, iview), t_chn(:), 'linear', 0);
            % add to image
            img_shot = img_shot + databuff;
        end
        % echo '.'
        if echo_onoff, fprintf('.'); end
    end
    % get img
    dataflow.image(Sxy ,imageindex) = gather(img_shot);
end

% normalize by Nviewprot, (pi/2 due to the filter)
dataflow.image = reshape(dataflow.image, imagesize, imagesize, Nimage).*(pi/Nviewprot/2);

% nolonger
% reorderflag = couchdirection < 0;
% dataflow.image = imagereorder(dataflow.image, Nslice, reorderflag);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function imageindex = getimageindex(Nslice, Nshot, ishot, couchdirection)

if couchdirection<0
    % backward couch
    imageindex = (1:Nslice) + (ishot-1)*Nslice;
else
    % forward couch
    imageindex = (1:Nslice) + (Nshot-ishot)*Nslice;
end

end