% % test data
% load F:\data-Dier.Z\3.2Head338\test\data_bnh_1122.mat

% prm
imagesize = single(prmflow.recon.imagesize);
h_img = single(prmflow.recon.FOV/prmflow.recon.imagesize);
delta_d = single(prmflow.recon.delta_d);
hond = h_img/delta_d;
Npixel = prmflow.recon.Npixel;
midchannel = prmflow.recon.midchannel;
Nviewprot = prmflow.recon.Nviewprot;
reconcenter = prmflow.recon.center;
delta_view = prmflow.recon.delta_view;
startviewangle = prmflow.recon.startviewangle;
Nimage = prmflow.recon.Nimage;
Nslice = single(prmflow.recon.Nslice);
FOV = prmflow.recon.FOV;
effFOV = single(min(FOV*1.2, prmflow.recon.maxFOV));
effNp = floor(effFOV/delta_d) + 1;
SID = prmflow.recon.SID;
Neighb = single(prmflow.recon.Neighb);
Nextslice = prmflow.recon.Nextslice;
Nedge = single(2);

% prepare
% view angle
% viewangle = mod((0:Nviewprot/2-1).*delta_view + startviewangle(1) + pi/2, pi*2);
viewangle = mod((0:Nviewprot-1).*delta_view + startviewangle(1) + pi/2, pi*2);
Nview = length(viewangle);
eta_C = reconcenter(1).*sin(viewangle) - reconcenter(2).*cos(viewangle);
indexstart = floor(midchannel + (-effFOV/2 + eta_C)./delta_d);
indexstart(indexstart<1) = 1;
ctrIdx = midchannel+eta_C./delta_d+1-indexstart;
channelpos = ((1:Npixel)'-midchannel).*delta_d;
channelindex = single(1:Npixel)';
maxR = effFOV/2/delta_d;

% XY
h = FOV/imagesize;
xygrid = single((-(imagesize-1)/2 : (imagesize-1)/2).*h);
[X, Y] = ndgrid(xygrid);
XY = [X(:) Y(:)] - reconcenter;

SID_h = SID/h_img;
% zz_samp = single(-Nslice+1:Nslice);

% Filter
filter = prmflow.recon.filter;
Hlen = length(filter);

% costheta = cos(viewangle);
% sintheta = sin(viewangle);

Nedge = 2;
Q1 = zeros(imagesize^2, Nslice+Nedge*2, 'single');
Q2 = zeros(imagesize^2, Nslice+Nedge*2, 'single');
index_img = repmat((1:imagesize^2)', 1, Nslice);

% test shot 1.5
ishot = 1;
img_index = (-Nslice/2-Nedge+1 : Nslice/2+Nedge) + Nslice*ishot;
% I know length img_index is Nextslice
image0 = dataflow.image(:,:,img_index);
image1 = zeros(imagesize, imagesize, Nslice, 'single');
image2 = zeros(imagesize, imagesize, Nslice, 'single');
% % test one view
% iview = 3;
for iview = 1:Nviewprot/2
    theta = viewangle(iview);
    % sintheta_iview = sintheta(iview);
    % costheta_iview = costheta(iview);
    
    % FP
    % interpXY
    channelpos2 = [channelpos; -channelpos]./h;
    [interpX, interpY, cs_view] = parallellinearinterp2D2(imagesize, imagesize, channelpos2, theta, reconcenter./h);
    interpY_rep = repmat(interpY, 1, 1, Nslice);
    interpX_rep = repmat(interpX, 1, 1, Nslice);
    % interpZ
    Eta = -interpX.*sin(theta) + interpY.*cos(theta);
    Zeta = interpX.*cos(theta) + interpY.*sin(theta);
    Zeta(Npixel+1:end, :) = -Zeta(Npixel+1:end, :);
    Zscale = sqrt(SID_h.^2 - Eta.^2);
    Zscale = (Zscale + Zeta)./Zscale;
    Zgrid = [1:Nslice/2 -Nslice/2+1:0];
    Zshift = [zeros(1, Nslice/2) ones(1, Nslice/2).*Nslice];
    interpZ = Zscale(:)*(Zgrid-1/2) + Zshift + Nedge + 1/2;
    interpZ = reshape(interpZ, Npixel*2, imagesize, Nslice);
    % project
    G = interp3(image0, interpY_rep, interpX_rep, interpZ, 'linear', 0);
    D0 = reshape(sum(G, 2).*(abs(cs_view)*h), Npixel, Nslice*2);
    
    % Filter
    Df = ifft(fft(D0, Hlen).*filter, 'symmetric');
    Df = Df(1:Npixel, :);
    % Df = reshape(Df, Npixel*2, Nslice);
    
    % BP
    % X-Y to Zeta-Eta
    Eta = -XY(:,1).*sin(theta) + XY(:,2).*cos(theta);
    Zeta = XY(:,1).*cos(theta) + XY(:,2).*sin(theta);
    
    t_chn1 = Eta./delta_d + midchannel;
    P1 = interp1(channelindex, Df(:, 1:2:end), t_chn1(:), 'linear', 0);
    t_chn2 = -Eta./delta_d + midchannel;
    P2 = interp1(channelindex, Df(:, 2:2:end), t_chn2(:), 'linear', 0);
    
    SetaD = sqrt(SID.^2 - Eta.^2);
    Zgrid = (1:Nslice)-1/2;
    Zscale_p = SetaD./(SetaD+Zeta);
    Zscale_n = SetaD./(SetaD-Zeta);
    Z1_p = Zscale_p(:)*Zgrid;
    Z1_n = Zscale_n(:)*Zgrid;
    s1_p = Z1_p<=Nslice/2;
    s1_n = Z1_n<=Nslice/2;
    Z2_p = Zscale_p(:)*(Zgrid-Nslice) + Nslice;
    Z2_n = Zscale_n(:)*(Zgrid-Nslice) + Nslice;
    
    Q1(:, Nedge+1:Nslice/2+Nedge) = P1(:, 1:Nslice/2);
    Q1(:, Nedge+Nslice/2+1:Nslice+Nedge) = P2(:, Nslice/2+1:Nslice);
    Q1(:, 1:2) = repmat(Q1(:, 3), 1, 2);
    Q1(:, end-1:end) = repmat(Q1(:, end-2), 1, 2);
    
    Q2(:, Nedge+1:Nslice/2+Nedge) = P2(:, 1:Nslice/2);
    Q2(:, Nedge+Nslice/2+1:Nslice+Nedge) = P1(:, Nslice/2+1:Nslice);
    Q2(:, 1:2) = repmat(Q2(:, 3), 1, 2);
    Q2(:, end-1:end) = repmat(Q2(:, end-2), 1, 2);
    
    interpZ1 = Z1_p.*s1_p + Z2_n.*~s1_p + Nedge + 1/2;
    interpZ2 = Z1_n.*s1_n + Z2_p.*~s1_n + Nedge + 1/2;
    
    image1 = image1 + reshape(interp2(Q1, interpZ1, index_img, 'linear'), imagesize, imagesize, Nslice);
    image2 = image2 + reshape(interp2(Q2, interpZ2, index_img, 'linear'), imagesize, imagesize, Nslice);
end