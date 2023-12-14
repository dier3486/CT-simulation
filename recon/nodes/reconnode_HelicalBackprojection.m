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
% Z-upsampling
Zupsamp = recon.Zupsamp;
Zupsampling = Zupsamp.Zupsampling;
% I know Zupsampling>=1
if Zupsampling>1
    ZupMatrix = Zupsamp.ZupMatrix;
else
    ZupMatrix = [];
end
Nslice_up = (Nslice-1)*Zupsampling+1;
ConeWeightScale = recon.ConeWeightScale;

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
Zend = Zview(Nview - Nviewskip);
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

% GPU Array
if GPUonoff
    [Cd, Nslice, XY, delta_d, midchannel, Zend, upsampgamma, filter, filtlen, Zupsampling, ZupMatrix, ConeWeightScale, Zscale] = ...
        putinGPU(Cd, Nslice, XY, delta_d, midchannel, Zend, upsampgamma, filter, filtlen, Zupsampling, ZupMatrix, ConeWeightScale, Zscale);
end

% prepare same values
sigma_z = (Nslice-1)/2/Cd/Zscale;
ConeWeightScale_Cz = Cd*ConeWeightScale*Zscale;

for iblk = 1:Nviewblock
    imageindex = recon.startimgbyblk(iblk):recon.endimgbyblk(iblk);
    Nimgperblk = recon.endimgbyblk(iblk) - recon.startimgbyblk(iblk) + 1;
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
        D = sqrt(1-Eta.^2);
        Phi = asin(Eta)./(pi*2);
        Zv = Zviewblk(iview);
        Zf = Zv - Phi.*Cd;
        
        % interp target on Z
        Tz = (Zgridblk-Zf)./(D+Zeta).*Zscale;
        % cone weight
        Wcone = ((Nslice-1)/2 - abs(Tz)).*ConeWeightScale;
        % while Wcone<0 set Wcone=0; while Wcone>1 set Wcone=1.
        Wcone = Wcone.*(Wcone>0) - (Wcone-1).*(Wcone>1);
        % shift Tz to the Z-upsampled position
        Tz = (Tz + (Nslice-1)/2).*Zupsampling + 1;
        
        % interp target on Eta
        Tchn = repmat(Eta./delta_d + midchannel, 1, Nimgperblk);
        
        % interpolation on the porjection field
        data_2 = interp2(datablk(:,:, iview), Tz, Tchn, 'linear', 0);
       
        % to normalize the cone weight
        % I know, 
        %   sigma_z = (Nslice-1)/2/Cd/Zscale; 
        %   ConeWeightScale_Cz = Cd*ConeWeightScale*Zscale;

        % a0 & Wa0, to consider the first related view on same side
        a0 = (Zgridblk-Zv)./Cd+Phi - (D+Zeta).*sigma_z;
        ceila0 = ceil(a0);
        Wa0 = 1 - (ceila0 - a0)./(D+Zeta).*ConeWeightScale_Cz;
        % I know Wa0 shall >=0
        Wa0 = Wa0.*(Wa0>0) + ceila0;
        Ca0 = ceil(-Zv/Cd);
        % while (ceila0 < Ca0) Wa0 shall = Ca0
        Wa0 = Wa0 + (Ca0-Wa0).*(ceila0 < Ca0);

        % b0 & Wb0, to consider the last related view on same side
        b0 = (Zgridblk-Zv)./Cd+Phi + (D+Zeta).*sigma_z;
        floorb0 = floor(b0);
        Wb0 = 1 - (b0 - floorb0)./(D+Zeta).*ConeWeightScale_Cz;
        % I know Wb0 shall >=0
        Wb0 = floorb0 - Wb0.*(Wb0>0);
        Cb0 = floor((Zend-Zv)/Cd);
        % while (floorb0 > Cb0) Wb0 shall = Cb0 + 1
        Wb0 = Wb0 + (Cb0 - Wb0).*(floorb0 > Cb0) + 1;

        % api & Wapi, to consider the first related view on opposite
        api = (Zgridblk-Zv)./Cd-Phi - (D-Zeta).*sigma_z - 1/2;
        ceilapi = ceil(api);
        Wapi = 1 - (ceilapi - api)./(D-Zeta).*ConeWeightScale_Cz;
        Wapi = Wapi.*(Wapi>0) + ceilapi;
        Capi = ceil(-Zv/Cd - 1/2);
        Wapi = Wapi + (Capi - Wapi).*(ceilapi < Capi);
        
        % bpi & Wbpi, to consider the last related view on opposite
        bpi = (Zgridblk-Zv)./Cd-Phi + (D-Zeta).*sigma_z - 1/2;
        floorbpi = floor(bpi);
        Wbpi = 1 - (bpi - floorbpi)./(D-Zeta).*ConeWeightScale_Cz;
        Wbpi = floorbpi - Wbpi.*(Wbpi>0);
        Cbpi = floor((Zend-Zv)/Cd - 1/2);
        Wbpi = Wbpi + (Cbpi - Wbpi).*(floorbpi > Cbpi) + 1;

        % normalization upon the above components
        Wcone = Wcone./(Wb0 - Wa0 + Wbpi - Wapi);

        % add to imageblk
        imageblk = imageblk + data_2.*Wcone;
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

