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

% BP parameters
nodeprm = prmflow.pipe.(status.nodename);
% % interp method (linear or 4points)
% if isfield(BPprm, 'interp')
%     interpmethod = BPprm.interp;
% else
%     interpmethod = '4points';
% end
% % not support yet
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
Npixel = recon.Npixel;
Nslice = recon.Nslice;
Nview = recon.Nview;
Nviewprot = recon.Nviewprot;
Nviewblock = recon.Nviewblock;
viewblock = recon.viewblock;
imagesize = recon.imagesize;
Nimage = recon.Nimage;
NactiveXY = recon.NactiveXY;
if isfield(recon, 'Nviewskip')
    Nviewskip = recon.Nviewskip;
else
    Nviewskip = 0;
end
Sxy =  recon.activeXY(:);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

% ini image
dataflow.image = zeros(imagesize^2, Nimage, 'single');

% prepare the varibles
Zview = single(((0:Nview-1)-Nviewskip).*(recon.pitchlength/Nviewprot/recon.delta_z));
% I know pitchlength = -couchspeed*rotationspeed, whare the rotationspeed is
% the secs per rotation.


Cp = single(recon.pitch);  % I know Cp = -couchspeed*rotationspeed/(Nslice*delta_z);
Cd = Cp*Nslice;
sigma_z = single((Nslice-1)/Nslice);
XY = recon.XY;
delta_d = recon.delta_d/recon.SID;
midchannel = recon.midchannel;
upsampling = recon.upsampling;
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
index_pixel = single(1:NactiveXY)';
Zend = Zview(Nview - Nviewskip);
%     Zspd = single(recon.delta_z/recon.pitchlength);
channelindex = single(1:Npixel_up)';
upsampgamma = recon.upsampgamma;
if isFilter
    filter = recon.filter;
    filtlen = length(filter);
else
    filter = [];
    filtlen = [];
end

% GPU Array
if GPUonoff
    [Cp, Cd, sigma_z, XY, delta_d, midchannel, channelindex, index_pixel, Zend, upsampgamma, filter, filtlen] = ...
        putinGPU(Cp, Cd, sigma_z, XY, delta_d, midchannel, channelindex, index_pixel, Zend, upsampgamma, filter, filtlen);
end

for iblk = 1:Nviewblock
    imageindex = recon.startimgbyblk(iblk):recon.endimgbyblk(iblk);
    Nimgperblk = recon.endimgbyblk(iblk) - recon.startimgbyblk(iblk) + 1;
    if iblk<Nviewblock
        viewindex = (1:viewblock) + (iblk-1)*viewblock + Nviewskip;
        Nviewperblk = gpuArray(viewblock);
    else  % iblk == Nviewblock
        viewindex = ((iblk-1)*viewblock+1 : Nview-Nviewskip*2) + Nviewskip;
        Nviewperblk = gpuArray(Nview - Nviewskip*2 - viewblock*(Nviewblock-1));
    end

    imageblk = zeros(NactiveXY, Nimgperblk, 'single', 'gpuArray');
    viewangle = gpuArray(recon.viewangle(viewindex));
%     viewindexblk = gpuArray(single((viewindex-1-Nviewskip)./Nviewprot));  % for pi-helical recon
    datablk = gpuArray(dataflow.rawdata(:, :, viewindex));
    Zviewblk = gpuArray(Zview(viewindex));
    Zgridblk = gpuArray(recon.Zgrid(imageindex));
    
    % up sampling
    if upsampling
        datablk = doubleup(reshape(datablk, Npixel, []), upsampgamma);
    end

    % filter
    if isFilter
        datablk_f = fft(reshape(datablk, Npixel_up, []), filtlen);
        datablk_f = ifft(datablk_f.*filter, 'symmetric');
        datablk = reshape(datablk_f(1:Npixel_up, :), Npixel_up, Nslice, Nviewperblk);
    end

    costheta = cos(viewangle);
    sintheta = sin(viewangle);
    % loop the views
    for iview = 1:Nviewperblk
        % X-Y to Zeta-Eta
        Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
        Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);
        % interp on channel direction
        t_chn = Eta./delta_d + midchannel;
        data_1 = interp1(channelindex, datablk(:, :, iview), t_chn(:), 'linear', 0);

        % Zeta-Eta to Z
        D = sqrt(1-Eta.^2);
        Phi = asin(Eta)./(pi*2);
        Zv = Zviewblk(iview);
        Zf = Zv - Phi.*Cd;
%         Zfpi = Zv + Phi.*Cd + Cd/2;   % for pi-helical recon
        Tz = (Zgridblk-Zf)./(D+Zeta);
        Sz = abs(Tz) <= (Nslice-1)/2;
        Tz = (Tz + (Nslice+1)/2).*Sz;
        % interp on Z
        data_2 = interp2(data_1, Tz, repmat(index_pixel, 1, Nimgperblk), 'linear', 0);

%         % pi-helical recon
%         ZetaonD = Zeta./(Zeta+D);
%         PiC = ((Zgridblk.*Zspd-viewindexblk(iview)).*(1-ZetaonD) - Phi.*ZetaonD).*4;
%         Spi = PiC>=-1 & PiC<1;
%         imageblk = imageblk + data_2.*Spi;
        
        % view numbers in cone
        a0 = max((Zgridblk-Zv)./Cd+Phi - (D+Zeta)./2./Cp.*sigma_z, -Zv/Cd);
        b0 = min((Zgridblk-Zv)./Cd+Phi + (D+Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd);
        n0 = floor(b0) - ceil(a0) + 1;
        api = max((Zgridblk-Zv)./Cd-Phi - (D-Zeta)./2./Cp.*sigma_z, -Zv/Cd) - 1/2;
        bpi = min((Zgridblk-Zv)./Cd-Phi + (D-Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd) - 1/2;
        npi = floor(bpi) - ceil(api) + 1;
        
        % no cone-weight yet

        % add to imageblk
        imageblk = imageblk + data_2./(n0+npi);
    end
    % add to dataflow.image
    dataflow.image(Sxy, imageindex) = dataflow.image(Sxy, imageindex) + gather(imageblk);
    % the reliable images' index are 1 : recon.startimgbyblk(iblk+1)-1

    % echo '.'
    if echo_onoff, fprintf('.'); end
end

% normalize and reshape
dataflow.image = reshape(dataflow.image, imagesize, imagesize, Nimage).*(pi/Nviewprot);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

