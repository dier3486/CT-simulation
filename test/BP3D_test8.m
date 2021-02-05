% BP test code
% for 3D Axial slope rebin, shots fill up
% half shot

% load('E:\data\simulation\TM\test\tilt_test3.mat');

% inputs are dataflow, prmflow
if exist('df0', 'var')
    dataflow = df0;
    clear df0;
end

if exist('pf0', 'var')
    prmflow = pf0;
    clear pf0;
end

% detector = prmflow.system.detector;
% SID = detector.SID;
% SDD = detector.SDD;
% delta_z = detector.hz_ISO;
SID = prmflow.recon.SID;
delta_z = prmflow.recon.delta_z;
Nshot = prmflow.recon.Nshot;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
startviewangle = prmflow.recon.startviewangle + pi/2;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d/SID;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nimage = prmflow.recon.Nimage;
% imageincrement = prmflow.recon.imageincrement;
% imageincrement = delta_z;
imageincrement = 0.6;   % imageincrement<delta_z

imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
gantrytilt = single(prmflow.recon.gantrytilt);

FOV = 300;
h = FOV/imagesize;
Rfov = FOV/2*(sqrt(2)+1)/2;
Nedge = floor((Nslice*delta_z/2 - (sqrt(SID^2-Rfov^2) - Rfov)/SID*(Nslice-1)/2*delta_z)/imageincrement) + 2;
Nedge = min(Nedge, Nslice/2);
Nextslice = Nslice + Nedge*2;
% Nextslice = Nslice*2;
gpuDevice;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);


% ini image
img = zeros(imagesize, imagesize, Nimage, 'single');

xygrid = gpuArray(single((-(imagesize-1)/2 : (imagesize-1)/2).*(h/SID)));
[X, Y] = ndgrid(xygrid);
X = X(:);
Y = Y(:);

% zgrid = gpuArray(single((-(Nslice-1)/2 : (Nslice-1)/2)));
% zgrid = gpuArray(single(1/2 : Nslice-1/2));
% Zgrid = gpuArray(single(-Nextslice/2+1/2 : Nextslice/2-1/2));

% slicegrid = (-(Nslice-1)/2 : (Nslice-1)/2).*delta_z;
% midslice = (Nslice+1)/2;
index_np = gpuArray([Nslice/2+1:Nslice  1:Nslice/2]);
% index_np = gpuArray([(Nslice+1:Nextslice)  repmat(Nextslice, 1, Nslice-Nextslice/2)  ...
%             ones(1,Nslice-Nextslice/2)  1:Nextslice-Nslice]);
% index_self = gpuArray([false(1, Nedge) true(1, Nslice) false(1, Nedge)]);
% index_self = gpuArray(Nedge+1:Nedge+Nslice);
% index_neib = gpuArray([1:Nedge  Nedge+Nslice+1:Nedge*2+Nslice]);
% Zgrid_self = gpuArray(single(-(Nslice-1)/2 : (Nslice-1)/2));
% Zgrid_neib = gpuArray(single([(Nslice+1)/2:Nslice-1/2  -Nslice+1/2:-(Nslice+1)/2]));

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% edge expand 
Nfill0 = 4;

% Nextslice = gpuArray(Nextslice);

% ini GPU buffer
Nviewprot_gpu = gpuArray(single(Nviewprot));
Nslice_gpu = gpuArray(single(Nslice));
Nedge_gpu = gpuArray(single(Nedge));
Nextslice_gpu = gpuArray(single(Nextslice));
Nfill0_gpu = gpuArray(single(Nfill0));
midchannel_gpu = gpuArray(midchannel);
SID_gpu = gpuArray(SID);
imagesize_gpu = gpuArray(single(imagesize));
% delta_z_norm = gpuArray(delta_z/imageincrement/SID);
delta_d = gpuArray(delta_d);
Eta = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
Zeta = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% D = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_0 = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_pi = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% n_0 = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% n_pi = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% s_tmp = false(imagesize*imagesize, 1, 'gpuArray');
% gap_self = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% gap_neib = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% kg_self_0 = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
% kg_neib_0 = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
% kg_self_pi = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
% kg_neib_pi = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
% s_self_0 = false(imagesize*imagesize, Nslice, 'gpuArray');
% s_neib_0 = false(imagesize*imagesize, Nslice, 'gpuArray');
% s_self_pi = false(imagesize*imagesize, Nslice, 'gpuArray');
% s_neib_pi = false(imagesize*imagesize, Nslice, 'gpuArray');
% s_gap_self = false(imagesize*imagesize, Nslice, 'gpuArray');
% s_gap_neib = false(imagesize*imagesize, Nslice, 'gpuArray');

% gap = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_0 = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_prev = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% s_0 = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_pi = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_gap = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_gapprev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapnext = false(imagesize, imagesize, Nslice, 'gpuArray');
Tz = zeros(imagesize*imagesize, Nextslice, 'single', 'gpuArray');
% Tz_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% h_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_chn = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_chninv = zeros(imagesize, imagesize, 'single', 'gpuArray');
data_0 = zeros(imagesize*imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
% data_0 = zeros(imagesize, imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
data_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');
% t_chn_floor = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_chn_alpha = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_chninv_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_z_floor = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% t_z_alpha = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
alpha = zeros(imagesize*imagesize*Nextslice, 1, 'single', 'gpuArray');
beta = zeros(imagesize*imagesize*Nextslice, 3, 'single', 'gpuArray');
gamma = zeros(imagesize*imagesize, Nextslice, 'single', 'gpuArray');
% t_z_index = zeros(imagesize^2*Nslice, 1, 'double', 'gpuArray');
% t_z_coeff = zeros(imagesize^2*Nslice, 1, 'single', 'gpuArray');

index_img = repmat((1:imagesize_gpu*imagesize_gpu)', 1, Nextslice);
img_shot = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% img_self = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
% img_neib = zeros(imagesize*imagesize, Nedge*2, 'single', 'gpuArray');

% Nshift1 = repelem(cat(3, Nslice_gpu/2-1, -Nslice_gpu/2), imagesize, imagesize, Nslice_gpu);
% Nshift1 = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
% Nshift2 = repelem([Nslice_gpu/2-1 -Nslice_gpu/2], 1, Nslice_gpu/2);
% Nshift3 = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
% Nshift3 = reshape(repelem([2 Nslice_gpu+Nfill0+2], 1, Nextslice/2), 1, 1, []); 

channelindex = gpuArray(single(1:Npixel)');



% z interp prepare
Nzs = 512;
Rs = min(500, FOV*sqrt(2))/SID;
zeta_samp = linspace(-Rs/2, Rs/2, Nzs);
eta_samp = linspace(-Rs/2, Rs/2, Nzs);
% [zeta_s, eta_s] = ndgrid(zeta_samp, eta_samp);
[zeta_s, eta_s] = meshgrid(zeta_samp, eta_samp);
t_samp1 = axialconehomeomorph(eta_s, zeta_s, Nslice, gantrytilt);
Nleft = -Nslice/2+1-Nfill0+3/2;
t_samp1(t_samp1<Nleft) = Nleft;
Nright = Nslice/2+Nfill0-3/2;
t_samp1(t_samp1>Nright) = Nright;
t_samp1 = t_samp1 + Nslice/2+Nfill0;

% % z interp prepare2
% Nrs = 512;
% Nts = 1024;
% r_samp = linspace(0, Rs/2, Nrs);
% theta_samp = single(linspace(0, pi*2, Nts));
% zeta_s = cos(theta_samp(:))*r_samp;
% eta_s = sin(theta_samp(:))*r_samp;
% t_samp2 = axialconehomeomorph(eta_s, zeta_s, Nslice, gantrytilt);

% to gpu
zeta_samp = gpuArray(zeta_samp);
eta_samp = gpuArray(eta_samp);
zz_samp = -Nslice_gpu+1:Nslice_gpu;
t_samp1 = gpuArray(t_samp1);
% r_samp = gpuArray(r_samp);
% theta_samp = gpuArray(theta_samp);
% t_samp2 = gpuArray(t_samp2);
z_target = repmat(-Nextslice_gpu/2+1 : Nextslice_gpu/2, imagesize*imagesize, 1);

Rxy = sqrt(X.^2+Y.^2);
thetaxy = atan2(Y, X);

% interp prepare
Nintp = 512;
alpha_intp = gpuArray(single(linspace(0, 1, Nintp+1)'));
index_intp = gpuArray(single(1:Nintp+1));
Nintp_gpu = gpuArray(single(Nintp));

% coeff
gamma_coeff1 = 0.6;
gamma_coeff2 = 1.4;
gamma_coeff1 = gpuArray(gamma_coeff1);
gamma_coeff2 = gpuArray(gamma_coeff2);

beta_intp = 1/2-sqrt(1+alpha_intp.*(1-alpha_intp).*4)./2;
t_intp = [(1+alpha_intp-beta_intp)./2  (alpha_intp+beta_intp)./2  ...
          (gamma_coeff1/4)./sqrt(1-alpha_intp.*(1-alpha_intp).*gamma_coeff2)];

% DX_intp = (1:imagesize_gpu^2)';
DZ_intp = gpuArray(single(1:(Nslice+Nfill0*2))');
DZ_intp2 = gpuArray(single(1:(Nslice/2+Nfill0))');

n20 = gpuArray(single(20));

% Nshot = 1;
for ishot = gpuArray(1:Nshot)
% for ishot = 1:1
tic;
%     if ishot==0
%         % first
%         imageindex = 1:Nslice/2;
%         viewangle = viewangle_prot + startviewangle(1);
%         viewshift = 0;
%     elseif ishot==Nshot
%         % last
%         imageindex = (1:Nslice/2) + (Nshot*Nslice-Nslice/2);
%         viewangle = viewangle_prot + startviewangle(Nshot);
%         viewshift = 0;
%     else
%         imageindex = (1:Nslice) + (ishot*Nslice-Nslice/2);
%         viewangle = viewangle_prot + startviewangle(ishot);
%         viewshift = round((startviewangle(ishot) - startviewangle(ishot+1))/(pi*2)*Nviewprot);
%     end

    imageindex = (1:Nextslice) + ((ishot-1)*Nslice-(Nextslice-Nslice)/2);
    gatherindex = 1:Nextslice;
    if ishot==1
        imageindex(1:(Nextslice-Nslice)/2) = [];
        gatherindex(1:(Nextslice-Nslice)/2) = [];
    end
    if ishot==Nshot
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        gatherindex(end-(Nextslice-Nslice)/2+1:end) = [];
    end
       

    viewangle = viewangle_prot + startviewangle(ishot);

    costheta = cos(viewangle);
    sintheta = sin(viewangle);
    
    % ini
    img_shot(:) = 0;
    % pi pair methed loop half rot
    for iview = 1:Nviewprot
%     for iview = 1:gpuArray(single(10))
%     for iview = 1:Nviewprot_gpu/2
        % .
%         pause(0.1);
        fprintf('.');
% tic;
        % X-Y to Eta-Zeta
        sintheta_iview = sintheta(iview);
        costheta_iview = costheta(iview);
% toc;
% tic;
        Eta = -X.*sintheta_iview + Y.*costheta_iview;
        Zeta = X.*costheta_iview + Y.*sintheta_iview;
        
%         Eta = -X.*sintheta_iview;
%         Zeta = X.*costheta_iview;
        % I know the +pi is -Eta and -Zeta
        % samples step on Z
%         Eta = Rxy.*sin(thetaxy - viewangle(iview));
%         Zeta = Rxy.*cos(thetaxy - viewangle(iview));
        
% toc;
% tic;
    % Tz
%     for ii = 1:Nextslice_gpu
%         Tz(:,ii) = interp2(zeta_samp, eta_samp, t_samp(:,:, ii+Nslice_gpu/2-Nedge_gpu+1), Zeta, Eta);
%     end
    Tz = interp3(zeta_samp, eta_samp, zz_samp, t_samp1, repmat(Zeta, 1, Nextslice_gpu), repmat(Eta, 1, Nextslice_gpu), z_target);
%     Tz = interp3(r_samp, theta_samp, zz_samp, t_samp2, repmat(Rxy(:), 1, Nextslice_gpu), repmat(theta_iview, 1, Nextslice_gpu), z_target);
% toc;
% tic;
    % get data 
    data_iview  = gpuArray(dataflow.rawdata(:, :, iview, ishot));
% toc;
% tic;
    % interp target (channel)
    t_chn = Eta./delta_d + midchannel_gpu;
    data_0(:, Nfill0_gpu+1:end-Nfill0_gpu) = interp1(channelindex, data_iview, t_chn(:), 'linear', 0);
% toc;
    
    % first shot and last shot
%     if ishot==1
%         data_0(:, 1:Nfill0) = repamt(data_0(:, Nfill0+1), 1, Nfill0);
%     end
%     if ishot==Nshot
%         data_0(:, end-Nfill0+1:end) = repamt(data_0(:, end-Nfill0), 1, Nfill0);
%     end

    data_0(:, 1:Nfill0_gpu) = repmat(data_0(: ,Nfill0_gpu+1), 1, Nfill0_gpu).*(ishot == 1);
    data_0(:, end-Nfill0_gpu+1:end) = repmat(data_0(:, end-Nfill0_gpu), 1, Nfill0_gpu).*(ishot == Nshot);
    
% tic;
    % interp on Z
    Tz_floor = floor(Tz(:));
    s_odd = mod(Tz_floor, 2);
    alpha = Tz(:) - Tz_floor;
    %         beta = sqrt(alpha+1.0);
    %         beta = interp1(alpha_intp, t_intp, alpha);
    beta = interp1(index_intp, t_intp, alpha.*Nintp_gpu+1);
    % I know beta = interp1(t_intp, alpha.*Nintp_gpu+1); but which lay in a matlab bug
    t_odd = reshape((Tz_floor+s_odd)./2 + beta(:,2).*s_odd + beta(:,1).*(1-s_odd), imagesize_gpu^2, []);
    t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), imagesize_gpu^2, []);
    data1 = interp2(data_0(:, 1:2:end), t_odd, index_img, 'linear', 0)./2 + ...
        interp2(data_0(:, 2:2:end), t_even, index_img, 'linear', 0)./2;
    data1 = data1 + interp2(conv2(data_0, [-1 2 -1], 'same'), Tz, index_img, 'linear', 0) ...
        .*reshape(beta(:,3), imagesize_gpu^2, []);
    
%     % test, linear interp
%     data1 = interp2(data_0, Tz, index_img, 'linear', 0);
% toc;
% tic;
    % add to image
    img_shot = img_shot + reshape(data1, imagesize_gpu, imagesize_gpu, []);
% toc;
%     pause(0.1);
    end
    % get img
    img(:,:,imageindex) = img(:,:,imageindex) + gather(img_shot(:,:,gatherindex));
    fprintf('\n');
toc;
end

img = img.*(pi/Nviewprot/2);

