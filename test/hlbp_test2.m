% run after reconnode_Helicalrebin
load F:\SINO\data\568\Bay4-568-20210326\PatientRawData\Chentao\Hel1.mat

% filter
prmflow.pipe.Filter.name = 'hann';
prmflow.pipe.Filter.freqscale = 1.5;
status.nodename = 'Filter';
[dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);


prmflow.protocol.reconFOV = 250;
prmflow.pipe.Backprojection.FOV = 250;
prmflow.pipe.Backprojection.viewblock = 576;
status.nodename = 'Backprojection';
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

% % start from here, these codes will be written in reconnode

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

% reshape rawdata
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

Zview = single(((0:Nview-1)-Nviewskip).*(recon.pitchlength/Nviewprot/recon.delta_z));
% I know pitchlength = -couchspeed*rotationspeed, whare the rotationspeed is
% the secs per rotation.


Cp = gpuArray(single(recon.pitch));  % I know Cp = -couchspeed*rotationspeed/(Nslice*delta_z);
Cd = Cp*Nslice;
sigma_z = gpuArray(single((Nslice-1)/Nslice));
XY = gpuArray(recon.XY);
delta_d = gpuArray(recon.delta_d/recon.SID);
midchannel = gpuArray(recon.midchannel);
channelindex = gpuArray(single(1:Npixel)');
index_pixel = gpuArray(single(1:NactiveXY)');

Zend = gpuArray(Zview(Nview - Nviewskip*2));
Zspd = gpuArray(single(recon.delta_z/recon.pitchlength));

% ini image
image = zeros(imagesize^2, Nimage, 'single');

for iblk = 1:Nviewblock
% for iblk = 1:6
    imageindex = recon.startimgbyblk(iblk):recon.endimgbyblk(iblk);
    Nimgperblk = recon.endimgbyblk(iblk) - recon.startimgbyblk(iblk) + 1;
    if iblk<Nviewblock
        viewindex = (1:viewblock) + (iblk-1)*viewblock + Nviewskip;
        Nviewperblk = gpuArray(viewblock);
    else  % iblk == Nviewblock
        viewindex = ((iblk-1)*viewblock+1 : Nview-Nviewskip) + Nviewskip;
        Nviewperblk = gpuArray(Nview - Nviewskip*2 - viewblock*(Nviewblock-1));
    end

    imageblk = zeros(NactiveXY, Nimgperblk, 'single', 'gpuArray');
    viewangle = gpuArray(recon.viewangle(viewindex));
    viewindexblk = gpuArray(single((viewindex-1-Nviewskip)./Nviewprot));
    datablk = gpuArray(dataflow.rawdata(:, :, viewindex));
    Zviewblk = gpuArray(Zview(viewindex));
    Zgridblk = gpuArray(recon.Zgrid(imageindex));

    costheta = cos(viewangle);
    sintheta = sin(viewangle);

    for iview = 1:Nviewperblk
        % X-Y to Zeta-Eta
        Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
        Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);
        % interp on channel direction
        t_chn = Eta./delta_d + midchannel;
        data_1 = interp1(channelindex, datablk(:, :, iview), t_chn(:), 'linear', 0);

        D = sqrt(1-Eta.^2);
        Phi = asin(Eta)./(pi*2);
        Zv = Zviewblk(iview);
        Zf = Zv - Phi.*Cd;
        Zfpi = Zv + Phi.*Cd + Cd/2;

        Tz = (Zgridblk-Zf)./(D+Zeta);
        Sz = abs(Tz) <= (Nslice-1)/2;
        Tz = (Tz + (Nslice+1)/2).*Sz;
        data_2 = interp2(data_1, Tz, repmat(index_pixel, 1, Nimgperblk), 'linear', 0);

%         ZetaonD = Zeta./(Zeta+D);
%         PiC = ((Zgridblk.*Zspd-viewindexblk(iview)).*(1-ZetaonD) - Phi.*ZetaonD).*4;
%         Spi = PiC>=-1 & PiC<1;
        
        a0 = max((Zgridblk-Zv)./Cd+Phi - (D+Zeta)./2./Cp.*sigma_z, -Zv/Cd);
        b0 = min((Zgridblk-Zv)./Cd+Phi + (D+Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd);
        n0 = floor(b0) - ceil(a0) + 1;
        %     n0 = floor((Zgrid(imgindex)-Zf)./Cd + (D+Zeta)./2./Cp) - ceil((Zgrid(imgindex)-Zf)./Cd - (D+Zeta)./2./Cp) + 1;

        %     bpi = min((Zgrid(imgindex)-Zfpi)./Cd + (D-Zeta)./2./Cp, (Zend-Zview(iview)-pi*Zspeed)/Cd);
        %     api = max((Zgrid(imgindex)-Zfpi)./Cd - (D-Zeta)./2./Cp, (-Zview(iview)-pi*Zspeed)/Cd);
        api = max((Zgridblk-Zv)./Cd-Phi - (D-Zeta)./2./Cp.*sigma_z, -Zv/Cd) - 1/2;
        bpi = min((Zgridblk-Zv)./Cd-Phi + (D-Zeta)./2./Cp.*sigma_z, (Zend - Zv)/Cd) - 1/2;
        npi = floor(bpi) - ceil(api) + 1;
        % add to imageblk
%         imageblk = imageblk + data_2.*Spi;
        imageblk = imageblk + data_2./(n0+npi);
    end
    image(Sxy, imageindex) = image(Sxy, imageindex) + gather(imageblk);
end

image = reshape(image, imagesize, imagesize, Nimage).*(pi/Nviewprot);