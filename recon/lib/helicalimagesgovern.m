function [Nimage, Nviewskip, Nblock, startimgbyblk, endimgbyblk, startimgbyblk_pi, endimgbyblk_pi, imageZgrid, Zviewshift] = ...
    helicalimagesgovern(Reff, Nview, Nviewprot, Nslice, Pitch, imagesperpitch, viewblock, Nimage, startview2image)
% A helical prepare sub function, which will find out the relations between the views with the images and how many images to
% reconstruct.
% [recon.Nimage, recon.Nviewskip, recon.Nviewblock, recon.startimgbyblk, recon.endimgbyblk, recon.Zgrid, recon.Zviewshift] = ...
%    helicalimagesgovern(Reff, recon.Nview, recon.Nviewprot, recon.Nslice, recon.pitch, recon.imagesperpitch, ...
%    recon.viewblock, Nimage, startview2image);
% The Nimage and startview2image can be skipped or set to [] in automatically getting.
% Reff is the effective radius of recon FOV (on SID), suggest: 
%   Reff = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2)/recon.SID;

% Nimage is preset images number
if nargin<8
    Nimage = []; 
end
% startview2image is the Z distance from the focal position of first view to the first image center (on recon.imageincrement)
% along the images increment direction
if nargin<9
    startview2image = [];
end
% NOTE: in helical recon we will use the recon.imageincrement as the unit of the distance in Z direction. Therefore everything
% will be normalized by imageincrement, e.g. we will use the imagesperpitch to present the pitch or pitchlength.
% And we will avoid to use the recon.pitch and recon.nominalpitch in underlying codes due to their possible ambiguity, as you
% read here as you know.

% the governing of the images related with the views
Nviewskip = floor(asin(Reff)*Nviewprot/pi/2);
Nvieweff = Nview - Nviewskip*2;

if isempty(startview2image)
    Zviewshift = 0;
else
    Zviewshift = -startview2image + Nviewskip*imagesperpitch/Nviewprot;
end

% I know, recon.pitch = -couchspeed*rotationspeed/(Nslice*delta_z);
% and recon.nominalpitch = abs(couchspeed*rotationspeed/(Nslice*imageincrement));
absPitch = abs(Pitch);

% Suggest: imagesperpitch = pitchlength/imageincrement; mostly it is nominalpitch*Nslice.

% pitchlength = Cp*Nslice*delta_z;
% imagesperpitch = pitchlength/imageincrement;    % Cd is the images' number per pitch

% all view recon condition
sigma_z = (Nslice-1)/Nslice;
phi0 = fzero(@(x) (sin(x)+cos(x)*sin(x)/sqrt(Reff^2-sin(x)^2)).*sigma_z-absPitch/pi, 0);
Z0 = phi0/(pi*2)*absPitch + cos(phi0)/2.*sigma_z + sqrt(Reff^2-sin(phi0)^2)/2.*sigma_z;
Next_0 = Z0/absPitch*Nviewprot;

% pi-line recon condition
phi_pi = fzero(@(x) (sin(x)*(x/pi - 1/2))/(Reff^2 - sin(x)^2)^(1/2) - (Reff^2 - sin(x)^2)^(1/2)/(pi*cos(x)) - ...
    (sin(x)*(Reff^2 - sin(x)^2)^(1/2)*(x/pi - 1/2))/cos(x)^2, 0);
Zpi = (sqrt(Reff^2-sin(phi_pi)^2)/cos(phi_pi)*(1/4-phi_pi/pi/2)+1/4);
Next_pi = Zpi*Nviewprot;

if isempty(Nimage)
    Nimage_all = round(Nvieweff/Nviewprot*imagesperpitch);
    imageZgrid = 0 : (Nimage_all-1);
else
    imageZgrid = 0 : (Nimage-1);
end

% start/end view of each image in pi-line recon
Vstart_pi = ceil((imageZgrid-Zviewshift).*(Nviewprot/imagesperpitch) - Next_pi);
Vend_pi = floor((imageZgrid-Zviewshift).*(Nviewprot/imagesperpitch) + Next_pi);

% available images are which can be reconstructed by pi-line
imgavail = (Vstart_pi>0) & (Vend_pi <= Nvieweff);
if ~isempty(Nimage) && any(~imgavail)
    warning('The configured images range is too large, which will be reset!');
end
% set the available images
Nimage = sum(imgavail);
imageZgrid = imageZgrid(imgavail); 
% re-set the available start/end view
Vstart_pi = Vstart_pi(imgavail);
Vend_pi = Vend_pi(imgavail);

% start/end view of each image in full recon
Vstart_0 = ceil((imageZgrid-Zviewshift).*(Nviewprot/imagesperpitch) - Next_0);
Vend_0 = floor((imageZgrid-Zviewshift).*(Nviewprot/imagesperpitch) + Next_0);

% to calculate the start/end image of each view block by the start/end view of each image
% I mean, we know which the views are related with an image, use it we can find out the set of images related by a block of
% views. Those a set of images are from the startimgbyblk(ii) to the endimgbyblk(ii) for the ii-th view block;
Nblock = ceil(Nvieweff/viewblock);
startimgbyblk = nan(1, Nblock);
endimgbyblk = nan(1, Nblock);
startimgbyblk_pi = nan(1, Nblock);
endimgbyblk_pi = nan(1, Nblock);
for ii = 1:Nimage
    % end image
    % full
    Vstart_ii = max(1, ceil(Vstart_0(ii)/viewblock));
    Vend_ii = min(Nblock, ceil(Vend_0(ii)/viewblock));
    endimgbyblk(Vstart_ii:Vend_ii) = ii;

    % pi-line
    VstartPi_ii = max(1, ceil(Vstart_pi(ii)/viewblock));
    VendPi_ii = min(Nblock, ceil(Vend_pi(ii)/viewblock));
    endimgbyblk_pi(VstartPi_ii:VendPi_ii) = ii;
    
    % start image
    jj = Nimage+1-ii;
    % full
    Vstart_jj = max(1, ceil(Vstart_0(jj)/viewblock));
    Vend_jj = min(Nblock, ceil(Vend_0(jj)/viewblock));
    startimgbyblk(Vstart_jj:Vend_jj) = jj;

    % pi-line
    VstartPi_jj = max(1, ceil(Vstart_pi(jj)/viewblock));
    VendPi_jj = min(Nblock, ceil(Vend_pi(jj)/viewblock));
    startimgbyblk_pi(VstartPi_jj:VendPi_jj) = jj;
end

end