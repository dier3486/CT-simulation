function [dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status)
% recon node, BP prepare, set FOV, image size, center (XYZ), tilt ... for the images
% [dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

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
BPprm = prmflow.pipe.(status.nodename);

prmflow.recon = commonbpprepare(prmflow.recon, prmflow.protocol, prmflow.system, BPprm);

% recon method
if isfield(BPprm, 'method') && ~isempty(BPprm.method)
    recon_method = BPprm.method;
else
    % default BP method
    switch lower(prmflow.recon.scan)
        case 'axial'
            recon_method = '2D';
        case 'helical'
            recon_method = '';
        otherwise
            recon_method = '';
    end
end
if ~strncmpi(recon_method, prmflow.recon.scan, length(prmflow.recon.scan))
    prmflow.recon.method = [prmflow.recon.scan recon_method];
else
    prmflow.recon.method = recon_method;
end
% I know a default recon method is Axial2D

% switch axial or helical
switch lower(prmflow.recon.scan)
    case 'axial'
        % image center and image number
        prmflow.recon.Nimage = prmflow.recon.Nslice * prmflow.recon.Nshot;
        [imagecenter, reconcenter_2DBP] = imagescenterintilt(prmflow.recon.center, prmflow.recon);
        prmflow.recon.imagecenter = imagecenter;
        prmflow.recon.reconcenter_2DBP = reconcenter_2DBP;
    case 'helical'
        % pitch
        prmflow.recon.pitchlength = -prmflow.protocol.couchspeed*prmflow.protocol.rotationspeed;
        prmflow.recon.pitch = prmflow.recon.pitchlength/(prmflow.recon.Nslice*prmflow.recon.delta_z);
        if isfield(prmflow.protocol, 'pitchhelical')
            prmflow.recon.nominalpitch = prmflow.protocol.pitchhelical;
        else
            prmflow.recon.nominalpitch = abs(prmflow.recon.pitch);
        end
        % helical prepare
        prmflow.recon = helicalprepare(prmflow.recon, BPprm);
        
    otherwise
        warning('The %s recon is not availabe!', prmflow.recon.scan);
        % no topo
end

% prepare for 3D BP
switch lower(prmflow.recon.method)
    case 'axial3d'
        prmflow.recon = axial3Dprepare(prmflow.recon, BPprm);
    otherwise
        % do nothing
        1;
end

% check rebin
if isfield(prmflow, 'rebin')
    if ~isfield(prmflow.rebin, 'issloperebin') || ~prmflow.rebin.issloperebin
        if prmflow.protocol.gantrytilt~=0
            % It is a mistake!
            warning(['The reconstruction no longer support previous Axialrebin when gantry tilting! ' ...
                'Please replace the reconnode Axialrebin by the Sloperebin.']);
        end
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [imagecenter, reconcenter_2DBP] = imagescenterintilt(reconcenter, recon)
% Cin is recon center; Cout is the rotation center on images
% 'small-step' style for tilt shots 

imagecenter = repmat(-reconcenter(:)', recon.Nimage, 1);
% Y shfit
% Yshift = -(recon.imageincrement*tan(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Yshift = -(recon.imageincrement*sin(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Yshift = fliplr(Yshift);
end
Yshift = repmat(Yshift(:), recon.Nshot, 1);
reconcenter_2DBP = imagecenter;
reconcenter_2DBP(:, 2) = reconcenter_2DBP(:, 2) + Yshift;
% Z shift
% Zshift = (recon.imageincrement*sec(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Zshift = recon.imageincrement.*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Zshift = fliplr(Zshift);
end
Zshift = Zshift(:) - (0:recon.Nshot-1).*(recon.imageincrement.*recon.Nslice).*recon.couchdirection;
Zshift = Zshift(:) - recon.startcouch;
imagecenter = [imagecenter Zshift];
reconcenter_2DBP = [reconcenter_2DBP Zshift];

end


function recon = axial3Dprepare(recon, BPprm)
% more prepare works for 3D Axial

% recon range
reconD = sqrt(sum((recon.FOV/2+abs(recon.center)).^2))*2;
% Neighb and Nextslice
if isfield(BPprm, 'Neighb')
    Neighb = BPprm.Neighb;
else
    Rfov = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2);
    Neighb = floor((recon.Nslice*recon.delta_z/2 - (sqrt(recon.SID^2-Rfov^2) - Rfov)/recon.SID*(recon.Nslice-1) ...
             /2*recon.delta_z)/recon.imageincrement) + 2;
end
Neighb = min(Neighb, recon.Nslice/2);
recon.Neighb = Neighb;
recon.Nextslice = recon.Nslice + Neighb*2;

% Zinterp
% defualt coeff
if isfield(BPprm, 'Gamma')
    Gamma = BPprm.Gamma;
else
    Gamma = [0.6, 1.4];
end
if isfield(BPprm, 'Nzsample')
    Nzsample = BPprm.Nzsample;
else
    Nzsample = [512 256];
end
recon.Zinterp = omiga4table(Gamma, Nzsample, recon.maxFOV, reconD, recon.SID, recon.Nslice, recon.gantrytilt);

end

function recon = helicalprepare(recon, BPprm)
% more prepare works for 3D Axial

Nslice = recon.Nslice;
Nviewprot = recon.Nviewprot;
Nview = recon.Nview;
if isfield(BPprm, 'viewblock') && ~isempty(BPprm.viewblock)
    viewblock = BPprm.viewblock;
else
    viewblock = Nviewprot;
end
recon.viewblock = viewblock;

Rf = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2)/recon.SID;
nviewskip = floor(asin(Rf)/recon.delta_view);
Nvieweff = Nview - nviewskip*2;

% I know pitch = -couchspeed*rotationspeed/(Nslice*delta_z);
Cp = abs(recon.pitch);
Cd = Cp*Nslice;    % = pitchlength/delta_z
sigma_z = (Nslice-1)/Nslice;

phi0 = fzero(@(x) (sin(x)+cos(x)*sin(x)/sqrt(Rf^2-sin(x)^2)).*sigma_z-Cp/pi, 0);
Z0 = (phi0/(pi*2)*Cp + cos(phi0)/2.*sigma_z + sqrt(Rf^2-sin(phi0)^2)/2.*sigma_z)*Nslice;
Next_0 = Z0/Cd*Nviewprot;

phi_pi = fzero(@(x) (sin(x)*(x/pi - 1/2))/(Rf^2 - sin(x)^2)^(1/2) - (Rf^2 - sin(x)^2)^(1/2)/(pi*cos(x)) - ...
    (sin(x)*(Rf^2 - sin(x)^2)^(1/2)*(x/pi - 1/2))/cos(x)^2, 0);
Zpi = (sqrt(Rf^2-sin(phi_pi)^2)/cos(phi_pi)*(1/4-phi_pi/pi/2)+1/4);
Next_pi = Zpi*Nviewprot;

Nimage_a = round(Nvieweff/Nviewprot*Cd);
index_imga = 1:Nimage_a;
Vstart_pi = ceil((index_imga-1).*(Nviewprot/Cd) - Next_pi);
Vend_pi = floor((index_imga-1).*(Nviewprot/Cd) + Next_pi);
imgavl = (Vstart_pi>0) & (Vend_pi <= Nvieweff);
% Vstart = Vstart(imgavl);
Nimage = sum(imgavl);

Vstart_0 = ceil((index_imga(imgavl)-1).*(Nviewprot/Cd) - Next_0);
Vend_0 = floor((index_imga(imgavl)-1).*(Nviewprot/Cd) + Next_0);

Nblock = ceil(Nvieweff/viewblock);
startimg = nan(1, Nblock);
endimg = nan(1, Nblock);
for ii = 1:Nimage
    Vstart_ii = max(1, ceil(Vstart_0(ii)/viewblock));
    Vend_ii = min(Nblock, ceil(Vend_0(ii)/viewblock));
    endimg(Vstart_ii:Vend_ii) = ii;
    
    jj = Nimage+1-ii;
    Vstart_jj = max(1, ceil(Vstart_0(jj)/viewblock));
    Vend_jj = min(Nblock, ceil(Vend_0(jj)/viewblock));
    startimg(Vstart_jj:Vend_jj) = jj;
end

% Nimage and available images by view blocks
recon.Nviewblock = Nblock;
recon.Nviewskip = nviewskip;
recon.Nimage = Nimage;
recon.startimgbyblk = startimg;
recon.endimgbyblk = endimg;

Zgrid = single(0:Nimage_a-1);
recon.Zgrid = Zgrid(imgavl);

end
