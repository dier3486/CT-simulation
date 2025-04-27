function [dataflow, prmflow, status] = reconnode_HelicalPiLineBP(dataflow, prmflow, status)
% recon node, Helical pi-line  BP (afrer reconnode_BPprepare)
% [dataflow, prmflow, status] = reconnode_HelicalPiLineBP(dataflow, prmflow, status);

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

isFilter = isfield(nodeprm, 'Filter');

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end
% echo onoff
echo_onoff = status.debug.echo_onoff;

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
% Z-upsampling
Zupsamp = recon.Zupsamp;
Zupsampling = Zupsamp.Zupsampling;
if Zupsampling>1
    ZupMatrix = Zupsamp.ZupMatrix;
else
    ZupMatrix = [];
end
Nslice_up = (Nslice-1)*Zupsampling+1;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, [], Nslice, Nview);
if recon.couchdirection>0
    % forward couch
    dataflow.rawdata = flip(dataflow.rawdata, 2);
end

% ini image
dataflow.image = zeros(imagesize(1)*imagesize(2), Nimage, 'single');

% prepare the parameters
% Cp = single(recon.pitch);  % I know Cp = -couchspeed*rotationspeed/(Nslice*delta_z);
Cd = recon.imagesperpitch;
Zscale = recon.imageincrement/recon.delta_z;
% sigma_z = single((Nslice-1)/Nslice)/Zscale;

Zview = single(((0:Nview-1)-Nviewskip).*(Cd/Nviewprot)) + recon.Zviewshift;
% Zview is the 'focal' Z-position of the views by image coordinate.
XY = recon.XY;

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

% is filtering?
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

% GPU Array
if GPUonoff
    [Cd, Nslice, XY, delta_d, midchannel, upsampgamma, filter, filtlen, Zupsampling, ZupMatrix, Zscale] = ...
        putinGPU(Cd, Nslice, XY, delta_d, midchannel, upsampgamma, filter, filtlen, Zupsampling, ZupMatrix, Zscale);
end


for iblk = 1:Nviewblock
    imageindex = recon.startimgbyblk_pi(iblk):recon.endimgbyblk_pi(iblk);
    Nimgperblk = recon.endimgbyblk_pi(iblk) - recon.startimgbyblk_pi(iblk) + 1;
%     imageindex = recon.startimgbyblk(iblk):recon.endimgbyblk(iblk);
%     Nimgperblk = recon.endimgbyblk(iblk) - recon.startimgbyblk(iblk) + 1;
    if iblk<Nviewblock
        viewindex = (1:viewblock) + (iblk-1)*viewblock + Nviewskip;
        Nviewperblk = viewblock;
    else  % iblk == Nviewblock
        viewindex = ((iblk-1)*viewblock+1 : Nview-Nviewskip*2) + Nviewskip;
        Nviewperblk = Nview - Nviewskip*2 - viewblock*(Nviewblock-1);
    end
    
    if GPUonoff
        imageblk = zeros(NactiveXY, Nimgperblk, 'single', 'gpuArray');
        viewangle = gpuArray(dataflow.rawhead.viewangle(viewindex));
        datablk = gpuArray(dataflow.rawdata(:, :, viewindex));
        Zviewblk = gpuArray(Zview(viewindex));
        Zgridblk = gpuArray(recon.Zgrid(imageindex));
    else
        imageblk = zeros(NactiveXY, Nimgperblk, 'single');
        viewangle = dataflow.rawhead.viewangle(viewindex);
        datablk = dataflow.rawdata(:, :, viewindex);
        Zviewblk = Zview(viewindex);
        Zgridblk = recon.Zgrid(imageindex);
    end
    
    % X upsampling
    if recon.upsampling
        datablk = doubleup(reshape(datablk, Npixel, []), upsampgamma);
    end

    % filter
    if isFilter
        datablk_f = fft(reshape(datablk, Npixel_up, []), filtlen);
        datablk_f = ifft(datablk_f.*filter, 'symmetric');
        datablk = reshape(datablk_f(1:Npixel_up, :), Npixel_up, Nslice, Nviewperblk);
    end

    % Z upsampling
    if Zupsampling>1
        datablk = reshape(permute(datablk, [1 3 2]), Npixel_up*Nviewperblk, Nslice);
        datablk = datablk*ZupMatrix;
        datablk = permute(reshape(datablk, Npixel_up, Nviewperblk, Nslice_up), [1 3 2]);
    end
    % sin cos
    costheta = cos(viewangle);
    sintheta = sin(viewangle);

    % loop the views
    for iview = 1:Nviewperblk
        % X-Y to Zeta-Eta
        Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
        Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);

        % Zeta-Eta to Z
        Deta = sqrt(1-Eta.^2);
        Phi = asin(Eta)./(pi*2);
        Zv = Zviewblk(iview);
        Zf = Zv - Phi.*Cd;
        
        % interp target on Z
        Tz = (Zgridblk-Zf)./(Deta+Zeta);
        % in Pi?
        PiC = (Tz.*Deta./Cd - Phi).*4;
        Spi = PiC>=-1 & PiC<1;
        % Z scale
        Tz = Tz.*Zscale + (Nslice+1)/2;

        % extrap (for big pitch) 
        % Tz(Tz<1) = 1;  Tz(Tz>Nslice) = Nslice;
        Tz = Tz.*(Tz>=1 & Tz<=Nslice) + (Tz<1).*1.0 + (Tz>Nslice).*Nslice;
        
        % shift Tz to the Z-upsampled position
        Tz = (Tz - 1).*Zupsampling + 1;
        
        % interp target on Eta
        Tchn = repmat(Eta./delta_d + midchannel, 1, Nimgperblk);
        
        % interpolation on the porjection field
        data_2 = interp2(datablk(:,:, iview), Tz, Tchn, 'linear', 0);
        
        % add to imageblk
        imageblk = imageblk + data_2.*Spi;

    end
    % add to dataflow.image
    dataflow.image(Sxy, imageindex) = dataflow.image(Sxy, imageindex) + gather(imageblk);
    % the reliable images' index are 1 : recon.startimgbyblk(iblk+1)-1

    % echo '.'
    if echo_onoff, fprintf('.'); end
end

% normalize and reshape
dataflow.image = reshape(dataflow.image, imagesize(2), imagesize(1), Nimage).*(pi/Nviewprot);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

