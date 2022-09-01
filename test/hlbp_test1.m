% % % run after hlrb_test1.m
% % % load D:\matlab\data\simu\TM\Hbp_sph2.mat 
% load F:\SINO\matlab\Data\simu\TM\Hfbp_Shp1.mat 

% filter
prmflow.pipe.Filter.name = 'hann';
prmflow.pipe.Filter.freqscale = 1.5;
status.nodename = 'Filter';
[dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);

prmflow.protocol.reconFOV = 250;
prmflow.pipe.Backprojection.FOV = 250;
status.nodename = 'Backprojection';
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

% to be written in prepare
SID = prmflow.recon.SID;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
Nview = prmflow.recon.Nview;
% stard/end view
Ncutview = floor(asin(FOV/2/SID)*(Nviewprot/pi/2));
startview = Ncutview+1;
endview = Nview - Ncutview;
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
dataflow.rawdata = dataflow.rawdata(:, startview:endview);
prmflow.recon.startviewangle = prmflow.recon.startviewangle + Ncutview/Nviewprot*pi*2;
prmflow.recon.Nview = Nview - Ncutview*2;

% recon parameters
SID = prmflow.recon.SID;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
Nview = prmflow.recon.Nview;
startviewangle = prmflow.recon.startviewangle + pi/2;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d/SID;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
% Neighb = prmflow.recon.Neighb;
% Nextslice = prmflow.recon.Nextslice;
% Nimage = prmflow.recon.Nimage;
reconcenter = prmflow.recon.center;
couchdirection = prmflow.recon.couchdirection;

recon = prmflow.recon;

rotationspeed = prmflow.protocol.rotationspeed;
couchspeed = prmflow.protocol.couchspeed;


dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

Cd = -rotationspeed*couchspeed/recon.delta_z;
Cp = Cd/Nslice;  % pitch
sigma_z = (Nslice-1)/Nslice;
Rf = FOV/2/SID;

phi0 = fzero(@(x) (sin(x)+cos(x)*sin(x)/sqrt(Rf^2-sin(x)^2)).*sigma_z-Cp/pi, 0);
Z0 = (phi0/(pi*2)*Cp + cos(phi0)/2.*sigma_z + sqrt(Rf^2-sin(phi0)^2)/2.*sigma_z)*Nslice;
Next_0 = Z0/abs(rotationspeed*couchspeed/recon.delta_z)*Nviewprot;

phi_pi = fzero(@(x) (sin(x)*(x/pi - 1/2))/(Rf^2 - sin(x)^2)^(1/2) - (Rf^2 - sin(x)^2)^(1/2)/(pi*cos(x)) - ...
    (sin(x)*(Rf^2 - sin(x)^2)^(1/2)*(x/pi - 1/2))/cos(x)^2, 0);
Zpi = (sqrt(Rf^2-sin(phi_pi)^2)/cos(phi_pi)*(1/4-phi_pi/pi/2)+1/4)*Cd;
Next_pi = Zpi/abs(rotationspeed*couchspeed/recon.delta_z)*Nviewprot;

Nimage_a = round(abs(Nview/Nviewprot*Cd));
index_imga = 1:Nimage_a;
Vstart_pi = ceil((index_imga-1).*(Nviewprot/Cd) - Next_pi);
Vend_pi = floor((index_imga-1).*(Nviewprot/Cd) + Next_pi);
imgavl = (Vstart_pi>0) & (Vend_pi <= Nview);
% Vstart = Vstart(imgavl);
Nimage = sum(imgavl);

Vstart_0 = ceil((index_imga(imgavl)-1).*(Nviewprot/Cd) - Next_0);
Vend_0 = floor((index_imga(imgavl)-1).*(Nviewprot/Cd) + Next_0);

startimg = nan(1, Nview);
endimg = nan(1, Nview);
indexPi = 0:Next_0*2;
for ii = 1:Nimage
    Vstart_ii = max(1, Vstart_0(ii));
    Vend_ii = min(Nview, Vend_0(ii));
    endimg(Vstart_ii:Vend_ii) = ii;

    jj = Nimage+1-ii;
    Vstart_jj = max(1, Vstart_0(jj));
    Vend_jj = min(Nview, Vend_0(jj));
    startimg(Vstart_jj:Vend_jj) = jj;
end
viewavl = find(~isnan(startimg));
startimg = gpuArray(single(startimg));
endimg = gpuArray(single(endimg));

Nslice = gpuArray(Nslice);

dataflow.image = zeros(512, 512, Nimage, 'single');

% XY
h = FOV/imagesize/SID;
xygrid = single((-(imagesize-1)/2 : (imagesize-1)/2).*h);
[X, Y] = ndgrid(xygrid);
XY = [X(:) Y(:)] - reconcenter./SID;

index_img = gpuArray(single(1:imagesize^2)');

Zgrid = single(0:Nimage_a-1);
Zgrid = Zgrid(imgavl);
Zgrid = gpuArray(Zgrid);



Zview = gpuArray(single((0:Nview-1).*(-rotationspeed*couchspeed/Nviewprot/recon.delta_z)));
Zspeed = gpuArray(single(-rotationspeed*couchspeed/pi/2/recon.delta_z));
Zend = Zview(end);


delta_d = gpuArray(delta_d);
midchannel = gpuArray(midchannel);
channelindex = single(1:Npixel)';
XY = gpuArray(XY);
viewangle = gpuArray(single((0:Nview-1).*(pi*2)./Nviewprot));

costheta = cos(viewangle + startviewangle);
sintheta = sin(viewangle + startviewangle);

image_out = zeros(512, 512, Nimage, 'single', 'gpuArray');
image2 = image_out;

% sliceindex = [1 1:Nslice Nslice];
sliceindex = 1:Nslice;

for iview = viewavl
% for iview = 100:Nviewprot:Nview
    imgindex = startimg(iview):endimg(iview);
    Nimg_iv = endimg(iview)-startimg(iview)+1;
%     imgindex = 1:Nimage;
%     Nimg_iv = Nimage;

    data_iview  = gpuArray(dataflow.rawdata(:, sliceindex, iview));
    sintheta_iview = sintheta(iview);
    costheta_iview = costheta(iview);

    % X-Y to Zeta-Eta
    Eta = -XY(:, 1).*sintheta_iview + XY(:, 2).*costheta_iview;
    Zeta = XY(:, 1).*costheta_iview + XY(:, 2).*sintheta_iview;
    % interp on channel direction
    t_chn = Eta./delta_d + midchannel;
    data_1 = interp1(channelindex, data_iview, t_chn(:), 'linear', 0);
    
    D = sqrt(1-Eta.^2);
    Phi = asin(Eta)./(pi*2);
    Zv = Zview(iview);
    Zf = Zv - Phi.*Cd;
    Zfpi = Zv + Phi.*Cd + Cd/2;

    Tz = (Zgrid(imgindex)-Zf)./(D+Zeta);
    Sz = abs(Tz) <= (Nslice-1)/2;
    Tz = (Tz + (Nslice+1)/2).*Sz;
    data_2 = interp2(data_1, Tz, repmat(index_img, 1, Nimg_iv), 'linear', 0);

    
    ZetaonD = Zeta./(Zeta+D);
    PiC = ((Zgrid(imgindex)./Zspeed-viewangle(iview)).*(1-ZetaonD) - asin(Eta).*ZetaonD).*(2/pi);
    Spi = PiC>=-1 & PiC<1;
    
%     b0 = min((Zgrid(imgindex)-Zf)./Cd + (D+Zeta)./2./Cp, (Zend - Zview(iview))/Cd);
%     a0 = max((Zgrid(imgindex)-Zf)./Cd - (D+Zeta)./2./Cp, -Zview(iview)/Cd);
    a0 = max((Zgrid(imgindex)-Zv)./Cd+Phi - (D+Zeta)./2./Cp.*sigma_z, -Zv/Cd);
    b0 = min((Zgrid(imgindex)-Zv)./Cd+Phi + (D+Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd);
    n0 = floor(b0) - ceil(a0) + 1;
%     n0 = floor((Zgrid(imgindex)-Zf)./Cd + (D+Zeta)./2./Cp) - ceil((Zgrid(imgindex)-Zf)./Cd - (D+Zeta)./2./Cp) + 1;

%     bpi = min((Zgrid(imgindex)-Zfpi)./Cd + (D-Zeta)./2./Cp, (Zend-Zview(iview)-pi*Zspeed)/Cd);
%     api = max((Zgrid(imgindex)-Zfpi)./Cd - (D-Zeta)./2./Cp, (-Zview(iview)-pi*Zspeed)/Cd);
    api = max((Zgrid(imgindex)-Zv)./Cd-Phi - (D-Zeta)./2./Cp.*sigma_z, -Zv/Cd) - 1/2;
    bpi = min((Zgrid(imgindex)-Zv)./Cd-Phi + (D-Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd) - 1/2;
    npi = floor(bpi) - ceil(api) + 1;
%     npi = floor((Zgrid(imgindex)-Zfpi)./Cd + (D-Zeta)./2./Cp) - ceil((Zgrid(imgindex)-Zfpi)./Cd - (D-Zeta)./2./Cp) + 1;

    image_out(:,:,imgindex) = image_out(:,:,imgindex) + reshape(data_2.*Spi, imagesize, imagesize, Nimg_iv);
    image2(:,:,imgindex) = image2(:,:,imgindex) + reshape(data_2./(n0+npi), imagesize, imagesize, Nimg_iv);
end

% image_out = image_out.*(pi/Nview_recon/2);
image_out = image_out.*(pi/Nviewprot);
image2 = image2.*(pi/Nviewprot);
