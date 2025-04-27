function recon = helicalprepare(recon, nodeprm)
% more prepare works for 3D Helical

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

% Z upsampling
recon.Zupsamp = Zupsamplingprepare(recon, nodeprm);

% Cone weight
if isfield(nodeprm, 'ConeWeightScale')
    recon.ConeWeightScale = nodeprm.ConeWeightScale;
else
    recon.ConeWeightScale = 1.0;
end

% imagesnumber per pitch
recon.imagesperpitch = recon.pitchlength/recon.imageincrement;

% the governing of the images related with the views
% effRadius is the reffective BP radius (to ISO) on SID
recon.effRadius = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2)/recon.SID;
if isfield(nodeprm, 'imagenumber') && ~isnan(nodeprm.imagenumber)
    recon.imagenumber = nodeprm.imagenumber;
else
    recon.imagenumber = [];
end
if isfield(nodeprm, 'startview2image')
    recon.startview2image = nodeprm.startview2image / recon.imageincrement;
else
    recon.startview2image = [];
end
% startview2image is the Z distance from the focal position of first view to the first image center along the images increment
% direction, whose unit is image number (/imageincrement).
if isfield(nodeprm, 'startview2bypass')
    recon.startview2bypass = nodeprm.startview2bypass / recon.imageincrement;
else
    recon.startview2bypass = [];
end
% startview2bypass is the Z distance from the focal position of first view to one if an image shall be located.

% views by the images
recon = helicalimagesgovern2(recon);

% imageStart (depends on Zviewshift)
recon.imageStart = (-recon.Nviewskip*recon.imagesperpitch/recon.Nviewprot + recon.Zviewshift)* ...
    recon.couchdirection;

% image center
Zshift = (recon.imageStart - (0:recon.Nimage-1).*recon.couchdirection).*recon.imageincrement - recon.startcouch;
recon.imagecenter = [repmat(-recon.center(:)', recon.Nimage, 1)  Zshift(:)];
% reconcenter is useless in helical 
% recon.reconcenter_2DBP = recon.imagecenter;

% forward projection 
SID_h = recon.SID/recon.voxelsize;
FPchannelpos_h = recon.forward.FPchannelpos./recon.voxelsize;
Phi = asin(FPchannelpos_h./SID_h)./(pi*2);
Zgrid_fp = -(recon.Nslice-1)/2 : (recon.Nslice-1)/2;
recon.forward.ZFPa = (sqrt(SID_h^2 - FPchannelpos_h.^2)./SID_h)*Zgrid_fp - Phi.*recon.imagesperpitch;
recon.forward.ZFPb = Zgrid_fp./SID_h;
% to use the parallellinearinterp3D3.m

end


function recon = helicalimagesgovern2(recon)
% A helical prepare sub function, which will find out the relations between the views with the images and how many images to
% reconstruct.

% NOTE: in helical recon we will use the recon.imageincrement as the unit of the distance in Z direction. Therefore everything
% will be normalized by imageincrement, e.g. we will use the imagesperpitch to present the pitch or pitchlength.
% And we will avoid to use the recon.pitch and recon.nominalpitch in underlying codes due to the possible ambiguities.

% I know, recon.pitch = -couchspeed*rotationspeed/(Nslice*delta_z);
% and recon.nominalpitch = abs(couchspeed*rotationspeed/(Nslice*imageincrement));
% Suggest: imagesperpitch = pitchlength/imageincrement; mostly it is nominalpitch*Nslice.

% all view recon condition
% sigma_z = (recon.Nslice-1)/recon.Nslice;
% Reff = recon.effRadius;
% % phi0 = fzero(@(x) (sin(x)+cos(x)*sin(x)/sqrt(Reff^2-sin(x)^2)).*sigma_z-absPitch/pi, 0);
% % Z0 = phi0/(pi*2)*absPitch + cos(phi0)/2.*sigma_z + sqrt(Reff^2-sin(phi0)^2)/2.*sigma_z;
% % fixed by (more robust)
% t0 = fzero(@(x) (Reff + sqrt((1-Reff^2*x^2)/(1-x^2)) )*x*sigma_z - absPitch/pi, 0);
% Z0 = asin(Reff*t0) * absPitch/(pi*2) + (sqrt(1-Reff^2*t0^2) + Reff*sqrt(1-t0^2)) * sigma_z/2;
% recon.viewextra_full = Z0/absPitch*recon.Nviewprot;
% 
% % pi-line recon condition
% t_pi = fzero(@(x)  (Reff^2 * x * sqrt(1-x^2)/(1-Reff^2*x^2) - x/sqrt(1-x^2))*pi/2 - 1, 0);
% Zpi = (1+Reff*sqrt( (1-t_pi^2)/(1-Reff^2*t_pi^2) ) )/4 - asin(t_pi*Reff)/pi/2;
% recon.viewextra_pi = Zpi*recon.Nviewprot;

% viewextra in all view and pi-line recon condition
[recon.viewextra_full, recon.viewextra_pi] = helicalviewexta(recon.Nslice, recon.pitch, recon.effRadius, recon.Nviewprot);

% Those fzero codes shall be replaced by a table in products.

% the governing of the images related with the views
recon.Nviewskip = floor(asin(recon.effRadius)*recon.Nviewprot/pi/2);

% active view number
recon.Nviewact = recon.Nview - recon.Nviewskip + ceil(recon.viewextra_full);
recon.Nviewact = min(recon.Nviewact, recon.viewnumber - recon.Nviewskip*2);
% mostly the Nviewact is recon.viewnumber - recon.Nviewskip*2

Cp = double(recon.Nviewprot/recon.imagesperpitch);
Nimage_all = round(recon.Nviewact/Cp);
% Nimage_all shall be big enough
imageZgrid = 0 : (Nimage_all-1);

if ~isempty(recon.startview2image)
    % user set start image
    Zviewshift = -recon.startview2image + recon.Nviewskip/Cp;
elseif ~isempty(recon.startview2bypass)
    % user set the image shall bypass
    k = floor(recon.startview2bypass - (recon.Nviewskip + ceil(recon.viewextra_pi))/Cp);
    Zviewshift = k - recon.startview2bypass + recon.Nviewskip/Cp;
else
    % default start image
    Zviewshift = -ceil(recon.viewextra_pi)*recon.imagesperpitch/recon.Nviewprot;
    % which will make the viewstart_pi(1) = 1.
end

% start/end view of each image in pi-line recon
viewstart_pi = ceil((imageZgrid - Zviewshift).*Cp - recon.viewextra_pi);
viewend_pi = floor((imageZgrid - Zviewshift).*Cp + recon.viewextra_pi);

% available images are which can be reconstructed by pi-line
imgavail = (viewstart_pi>0) & (viewend_pi <= recon.Nviewact);
recon.Nimage = sum(imgavail);
if ~isempty(recon.imagenumber)
    if recon.imagenumber < recon.Nimage
        k = find(imgavail, 1) + recon.imagenumber;
        imgavail(k : end) = false;
        recon.Nimage = recon.imagenumber;
    elseif any(~imgavail)
        warning('The configured images range is too large, which will be reset!');
    end
else
    % set imagenumber
    recon.imagenumber = floor( floor((recon.viewnumber - double(recon.Nviewskip)*2) - double(recon.viewextra_pi)) / Cp ...
         + Zviewshift) + 1;
    % the recon.viewnumber and recon.imagenumber could be inf
end

% set the available images
imageZgrid = imageZgrid(imgavail); 
% re-set the available start/end view
viewstart_pi = viewstart_pi(imgavail);
viewend_pi = viewend_pi(imgavail);

% start/end view of each image in full recon
viewstart_full = ceil((imageZgrid - Zviewshift).*Cp - recon.viewextra_full);
viewend_full = floor((imageZgrid - Zviewshift).*Cp + recon.viewextra_full);

% % Zviewend
% recon.Zviewend = (recon.Nviewact-1)*(recon.imagesperpitch/recon.Nviewprot);
% % I know Zviewstart = 0
% the [Zviewstart Zviewend] was replaced by [0 recon.Nviewact-1]

% ZviewRange
recon.ZviewRange = [0 recon.Nviewact-1];

% imagerely
recon.imagerely = ceil(recon.viewextra_full / Cp);

% to return the start/end view and other
recon.viewbyimages_pi = [viewstart_pi; viewend_pi];
recon.viewbyimages_full = [viewstart_full; viewend_full];
recon.imageZgrid = imageZgrid;
recon.Zviewshift = Zviewshift;      
% We suggest ot use the Zviewshift by adding it on Zview, but not adding it on imageZgrid due to the imageZgrid are integers.

% update Nview (Nviewskip has been read) 
recon.viewnumber = recon.viewnumber - recon.Nviewskip;
recon.Nview = recon.Nview - recon.Nviewskip;

end