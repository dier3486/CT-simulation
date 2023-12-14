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
recon = prmflow.recon;
SID = recon.SID;
Nshot = recon.Nshot;
% FOV = recon.FOV;
Nviewprot = recon.Nviewprot;
% startviewangle = recon.startviewangle;
imagesize = recon.imagesize;
% midchannel = recon.midchannel;
% delta_d = recon.delta_d/SID;
Npixel = recon.Npixel;
Nslice = recon.Nslice;
Nimage = recon.Nimage;
% reconcenter = recon.center;
couchdirection = recon.couchdirection;
viewblock = recon.viewblock;
Nviewblock = ceil(Nviewprot/viewblock);
XY = recon.XY;
Sxy = recon.activeXY;
NactiveXY = recon.NactiveXY;
% viewangle = recon.viewangle;

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

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, [], Nslice, Nviewprot, Nshot);
if couchdirection>0
    % forward couch
    dataflow.rawdata = flip(dataflow.rawdata, 2);
end
% ini image
dataflow.image = zeros(imagesize(1)*imagesize(2), Nimage, 'single');
% channel interp
channelindex = single(1:Npixel_up)';

% to GPU
if GPUonoff
    [XY, channelindex, delta_d, midchannel, upsampgamma, filter, filtlen, Npixel, Npixel_up] = ...
        putinGPU(XY, channelindex, delta_d, midchannel, upsampgamma, filter, filtlen, Npixel, Npixel_up);
end

for ishot = 1 : Nshot
    imageindex = (1:Nslice) + (ishot-1)*Nslice;
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
        viewangleblk = viewangle(viewindex);
        % get data per block
        if GPUonoff
            datablk = gpuArray(dataflow.rawdata(:, :, viewindex, ishot));
            [viewangleblk, Nviewperblk]  = putinGPU(viewangleblk, Nviewperblk);
        else
            datablk = dataflow.rawdata(:, :, viewindex, ishot);
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
dataflow.image = reshape(dataflow.image, imagesize(2), imagesize(1), Nimage).*(pi/Nviewprot/2);

% nolonger
% reorderflag = couchdirection < 0;
% dataflow.image = imagereorder(dataflow.image, Nslice, reorderflag);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


% function imageindex = getimageindex(Nslice, Nshot, ishot, couchdirection)
% 
% if couchdirection<0
%     % backward couch
%     imageindex = (1:Nslice) + (ishot-1)*Nslice;
% else
%     % forward couch
%     imageindex = (1:Nslice) + (Nshot-ishot)*Nslice;
% end
% 
% end