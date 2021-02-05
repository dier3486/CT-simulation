% BP test code
% for 3D Axial slope rebin, shots fill up
% half shot

% load('E:\data\simulation\TM\test\vol_test1.mat');

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
delta_d = prmflow.recon.delta_d;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nimage = prmflow.recon.Nimage;
% imageincrement = prmflow.recon.imageincrement;
% imageincrement = delta_z;
imageincrement = 0.6;   % imageincrement<delta_z

imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
% set center to (0, 0)
% imagecenter(:, 1:2) = 0;
% no tilt

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

xygrid = gpuArray(single((-(imagesize-1)/2 : (imagesize-1)/2).*h));
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
Zgrid_self = gpuArray(single(-(Nslice-1)/2 : (Nslice-1)/2));
Zgrid_neib = gpuArray(single([(Nslice+1)/2:Nslice-1/2  -Nslice+1/2:-(Nslice+1)/2]));

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% edge expand 
Nfill0 = gpuArray(single(4));
% Nextslice = gpuArray(Nextslice);

% ini GPU buffer
Nviewprot_gpu = gpuArray(single(Nviewprot));
Nslice_gpu = gpuArray(single(Nslice));
Nedge_gpu = gpuArray(single(Nedge));
midchannel_gpu = gpuArray(midchannel);
SID_gpu = gpuArray(SID);
imagesize_gpu = gpuArray(single(imagesize));
delta_z_norm = gpuArray(delta_z/imageincrement/SID);
Eta = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
Zeta = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
D = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
t_0 = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
t_pi = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
n_0 = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
n_pi = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
s_tmp = false(imagesize*imagesize, 1, 'gpuArray');
gap_self = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
gap_neib = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
kg_self_0 = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
kg_neib_0 = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
kg_self_pi = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
kg_neib_pi = zeros(imagesize*imagesize, Nslice, 'single', 'gpuArray');
s_self_0 = false(imagesize*imagesize, Nslice, 'gpuArray');
s_neib_0 = false(imagesize*imagesize, Nslice, 'gpuArray');
s_self_pi = false(imagesize*imagesize, Nslice, 'gpuArray');
s_neib_pi = false(imagesize*imagesize, Nslice, 'gpuArray');
s_gap_self = false(imagesize*imagesize, Nslice, 'gpuArray');
s_gap_neib = false(imagesize*imagesize, Nslice, 'gpuArray');

% gap = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_0 = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_prev = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% s_0 = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_pi = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_gap = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_gapprev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapnext = false(imagesize, imagesize, Nslice, 'gpuArray');
% Tz_0 = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% Tz_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% h_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_chn = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv = zeros(imagesize, imagesize, 'single', 'gpuArray');
data_0 = zeros(imagesize*imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
% data_0 = zeros(imagesize, imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
data_iview = zeros(Npixel, Nslice, 2, 'single', 'gpuArray');
t_chn_floor = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
t_chn_alpha = zeros(imagesize*imagesize, 1, 'single', 'gpuArray');
% t_chninv_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_z_floor = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% t_z_alpha = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
alpha = zeros(imagesize*imagesize*Nslice_gpu, 1, 'single', 'gpuArray');
beta = zeros(imagesize*imagesize*Nslice_gpu, 3, 'single', 'gpuArray');
gamma = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
% t_z_index = zeros(imagesize^2*Nslice, 1, 'double', 'gpuArray');
% t_z_coeff = zeros(imagesize^2*Nslice, 1, 'single', 'gpuArray');

index_img = repmat((1:imagesize_gpu*imagesize_gpu)', 1, Nslice_gpu);
% img_shot = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
img_self = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
img_neib = zeros(imagesize*imagesize, Nedge*2, 'single', 'gpuArray');

% Nshift1 = repelem(cat(3, Nslice_gpu/2-1, -Nslice_gpu/2), imagesize, imagesize, Nslice_gpu);
Nshift1 = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
Nshift2 = repelem([Nslice_gpu/2-1 -Nslice_gpu/2], 1, Nslice_gpu/2);
Nshift3 = zeros(imagesize*imagesize, Nslice_gpu, 'single', 'gpuArray');
% Nshift3 = reshape(repelem([2 Nslice_gpu+Nfill0+2], 1, Nextslice/2), 1, 1, []); 

channelindex = gpuArray(single(1:Npixel)');

% coeff
gamma_coeff1 = 0.6;
gamma_coeff2 = 1.5;
gamma_coeff1 = gpuArray(gamma_coeff1);
gamma_coeff2 = gpuArray(gamma_coeff2);

% interp prepare
Nintp = 512;
alpha_intp = gpuArray(single(linspace(0, 1, Nintp+1)'));
index_intp = gpuArray(single(1:Nintp+1));
Nintp_gpu = gpuArray(single(Nintp));

beta_intp = 1/2-sqrt(1+alpha_intp.*(1-alpha_intp).*4)./2;
t_intp = [(1+alpha_intp-beta_intp)./2  (alpha_intp+beta_intp)./2  ...
          gamma_coeff1./sqrt(1-alpha_intp.*(1-alpha_intp).*gamma_coeff2)];

% DX_intp = (1:imagesize_gpu^2)';
DZ_intp = gpuArray(single(1:(Nslice+Nfill0*2))');
DZ_intp2 = gpuArray(single(1:(Nslice/2+Nfill0))');

for ishot = gpuArray(1:Nshot)
% for ishot = 1:1
% tic;
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
    for iview = 1:Nviewprot/2
%     for iview = 1:Nviewprot_gpu/2
        % .
        fprintf('.');
tic;
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
toc;
tic;
        D = sqrt(SID_gpu^2 - Eta.^2);
        t_0 = (D + Zeta).*delta_z_norm;
        t_pi = (D - Zeta).*delta_z_norm;
%         gap = Nslice_gpu - (t_0+t_pi).*(Nslice_gpu-1)./2;
% toc;
% tic;        
        n_0 = floor((Nslice_gpu - t_pi.*(Nslice_gpu-1)./2)./t_0 - 0.5) + 1;
%         n_0(n_0>Nslice_gpu/2) = Nslice_gpu/2;
        s_tmp = n_0>Nslice_gpu/2;
        n_0 = n_0.*~s_tmp + (Nslice_gpu/2).*s_tmp;
        n_pi = floor((Nslice_gpu - t_0.*(Nslice_gpu-1)./2)./t_pi - 0.5) + 1;
        s_tmp = n_pi>Nslice_gpu/2;
        n_pi = n_pi.*~s_tmp + (Nslice_gpu/2).*s_tmp;
% toc;
% tic;
%         gap_0 = Nslice_gpu - t_0.*((Nslice_gpu-1)/2) - t_pi.*(n_pi-0.5);
%         gap_pi = Nslice_gpu - t_pi.*((Nslice_gpu-1)/2) - t_0.*(n_0-0.5);
        gap_self = Nslice_gpu - t_0.*((Nslice_gpu-1)/2) - t_pi.*(n_pi-0.5);
        gap_neib = Nslice_gpu - t_pi.*((Nslice_gpu-1)/2) - t_0.*(n_0-0.5);
        
        s_tmp = gap_self<eps;
        gap_self = gap_self + eps.*s_tmp;
        s_tmp = gap_neib<eps;
        gap_neib = gap_neib + eps.*s_tmp;
% toc;        
% tic;
        % samples from '0' and 'pi' scan
%         kg_0 = reshape(Zgrid,1,1,[])./t_0 + 1/2;
%         kg_pi = reshape(Zgrid,1,1,[])./t_pi + 1/2;
        
        kg_self_0 = Zgrid_self./t_0 + 1/2;
        kg_self_pi = Zgrid_self./t_pi + 1/2;
        kg_neib_0 = Zgrid_neib./t_0 + 1/2;
        kg_neib_pi = Zgrid_neib./t_pi + 1/2;
        
% toc;      
% tic;
%         s_0(:,:,index_self) = (kg_0(:,:,index_self) <= Nslice_gpu/2) & (kg_0(:,:,index_self) >= -Nslice_gpu/2+1);
%         s_0(:,:,~index_self) = (kg_0(:,:,~index_self) <= n_0) & (kg_0(:,:,~index_self) >= -n_0+1);
        s_self_0 = (kg_self_0 <= Nslice_gpu/2) & (kg_self_0 >= -Nslice_gpu/2+1);
        s_neib_0 = (kg_neib_0 <= n_0) & (kg_neib_0 >= -n_0+1);
        
%         s_pi(:,:,index_self) = (kg_pi(:,:,index_self) <= Nslice_gpu/2) & (kg_pi(:,:,index_self) >= -Nslice_gpu/2+1);
%         s_pi(:,:,~index_self) = (kg_pi(:,:,~index_self) <= n_pi) & (kg_pi(:,:,~index_self) >= -n_pi+1);
        s_self_pi = (kg_self_pi <= Nslice_gpu/2) & (kg_self_pi >= -Nslice_gpu/2+1);
        s_neib_pi = (kg_neib_pi <= n_pi) & (kg_neib_pi >= -n_pi+1);

% %         s_0 = (kg_0 < Nslice_gpu/2) & (kg_0 > -Nslice_gpu/2+1);   
% %         s_pi = (kg_pi < Nslice_gpu/2) & (kg_pi > -Nslice_gpu/2+1);
%         % I know s_neib_0 = s_self_0(index_np)
% toc;
% tic;
        s_gap_self = ~s_self_0 & ~s_neib_pi;
        s_gap_neib = ~s_neib_0 & ~s_self_pi;
%         s_gap = ~s_0 & ~s_pi(:,:,index_np);
%         % I know s_gap_pi = s_gap_0(index_np)
%         
% toc;
% tic;
%         % t_z
        Nshift1 = cat(2, repmat(-n_0-Nslice_gpu/2, 1, Nslice_gpu/2), repmat(n_0+Nslice_gpu/2, 1, Nslice_gpu/2));
%         Nshift2 = 0;
        Tz_self_0 = kg_self_0.*s_self_0 + (kg_neib_pi+Nshift1).*s_neib_pi + ((kg_self_0 + Nshift2).*(t_0./gap_self) - Nshift2).*s_gap_self + (Nslice_gpu/2+Nfill0);
%         Nshift1 = cat(3, repmat(n_0+Nslice_gpu/2, 1, 1, Nslice_gpu/2), repmat(-n_0+Nslice_gpu/2, 1, 1, Nslice_gpu/2));
        Nshift3 = cat(2, repmat(-n_0, 1, Nslice_gpu/2), repmat(n_0-1, 1, Nslice_gpu/2));
        Tz_neib_0 = kg_neib_0.*s_neib_0 + (kg_self_pi+Nshift1(:,index_np)).*s_self_pi + ((kg_neib_0 + Nshift3).*(t_0./gap_neib) - Nshift3).*s_gap_neib + (Nslice_gpu/2+Nfill0);
        
        Nshift1 = cat(2, repmat(-n_pi-Nslice_gpu/2, 1, Nslice_gpu/2), repmat(n_pi+Nslice_gpu/2, 1, Nslice_gpu/2));
        Tz_self_pi = kg_self_pi.*s_self_pi + (kg_neib_0+Nshift1).*s_neib_0 + ((kg_self_pi + Nshift2).*(t_pi./gap_neib) - Nshift2).*s_gap_neib + (Nslice_gpu/2+Nfill0);
        Nshift3 = cat(2, repmat(-n_pi, 1, Nslice_gpu/2), repmat(n_pi-1, 1, Nslice_gpu/2));
        Tz_neib_pi = kg_neib_pi.*s_neib_pi + (kg_self_0+Nshift1(:,index_np)).*s_self_0 + ((kg_neib_pi + Nshift3).*(t_pi./gap_self) - Nshift3).*s_gap_self + (Nslice_gpu/2+Nfill0);
        
% %         Tz_0(:,:,1:Nslice_gpu) = -Nslice_gpu/2;     Tz_0(:,:,Nslice_gpu+1:end) = Nslice_gpu/2;
% %         Tz_pi(:,:,1:Nslice_gpu) = -Nslice_gpu/2;     Tz_pi(:,:,Nslice_gpu+1:end) = Nslice_gpu/2;
%         
% %         Nshift1 = cat(3, repmat(n_0-1, 1, 1, Nextslice/2), repmat(-n_0, 1, 1, Nextslice/2));
%         Nshift1(:,:,1:Nedge) = repmat(n_0-1, 1, 1, Nedge);  Nshift1(:,:,end-Nedge+1:end) = repmat(-n_0, 1, 1, Nedge);
%         Tz_0 = (kg_pi(:,:,index_np) + Nshift2).*s_pi(:,:,index_np) + kg_0.*s_0 + ((kg_0 + Nshift1).*(t_0./gap) - Nshift1).*s_gap + Nslice_gpu/2+Nfill0;
%         Nshift1(:,:,1:Nedge) = repmat(n_pi-1, 1, 1, Nedge);  Nshift1(:,:,end-Nedge+1:end) = repmat(-n_pi, 1, 1, Nedge);
%         Tz_pi = (kg_0(:,:,index_np) + Nshift2).*s_0(:,:,index_np) + kg_pi.*s_pi + ((kg_pi + Nshift1).*(t_pi./gap(:,:,index_np)) - Nshift1).*s_gap(:,:,index_np) + Nslice_gpu/2+Nfill0;
% toc;
% tic;
        
        s_tmp1 = Tz_self_0<2;
        s_tmp2 = Tz_self_0>Nslice_gpu+Nfill0+2;
        Tz_self_0 = Tz_self_0.*(~s_tmp1 & ~s_tmp2) + s_tmp1.*2 + s_tmp2.*(Nslice_gpu+Nfill0+2);
        
        s_tmp1 = Tz_self_pi<2;
        s_tmp2 = Tz_self_pi>Nslice_gpu+Nfill0+2;
        Tz_self_pi = Tz_self_pi.*(~s_tmp1 & ~s_tmp2) + s_tmp1.*2 + s_tmp2.*(Nslice_gpu+Nfill0+2);
        
        Tz_neib_0 = [Tz_neib_0(:, 1:Nedge_gpu) Tz_neib_0(:, end-Nedge_gpu+1:end)];
        s_tmp1 = Tz_neib_0<2;
        s_tmp2 = Tz_neib_0>Nslice_gpu+Nfill0+2;
        Tz_neib_0 = Tz_neib_0.*(~s_tmp1 & ~s_tmp2) + s_tmp1.*2 + s_tmp2.*(Nslice_gpu+Nfill0+2);
        
        Tz_neib_pi = [Tz_neib_pi(:, 1:Nedge_gpu) Tz_neib_pi(:, end-Nedge_gpu+1:end)];
        s_tmp1 = Tz_neib_pi<2;
        s_tmp2 = Tz_neib_pi>Nslice_gpu+Nfill0+2;
        Tz_neib_pi = Tz_neib_pi.*(~s_tmp1 & ~s_tmp2) + s_tmp1.*2 + s_tmp2.*(Nslice_gpu+Nfill0+2);
        
%         Tz_self_0(Tz_self_0<2) = 2;    Tz_self_0(Tz_self_0>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
%         Tz_self_pi(Tz_self_pi<2) = 2;    Tz_self_pi(Tz_self_pi>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
        
%         Tz_0(Tz_0<2) = 2;    Tz_0(Tz_0>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
%         Tz_pi(Tz_pi<2) = 2;    Tz_pi(Tz_pi>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
        
%         s_0 = (Tz_0 >= 2) & (Tz_0 <= (Nslice_gpu+Nfill0+2));
%         Tz_0 = Tz_0.*s_0 + Nshift3.*~s_0;
%         s_pi = (Tz_pi >= 2) & (Tz_pi <= (Nslice_gpu+Nfill0+2));
%         Tz_pi = Tz_pi.*s_pi + Nshift3.*~s_pi;
        
%         s_0 = (Tz_0 >= 2);
%         Tz_0(:,:,1:Nslice_gpu) = (Tz_0 - 2).*s_0 + 2;
%         s_0 = Tz_0 <= (Nslice_gpu+Nfill0+2);
%         Tz_0 = (Tz_0 - (Nslice_gpu+Nfill0+2)).*s_0 + (Nslice_gpu+Nfill0+2);
        
%         s_pi = Tz_pi >= 2;
%         Tz_pi = (Tz_pi - 2).*s_pi + 2;
%         s_pi = Tz_pi <= (Nslice_gpu+Nfill0+2);
%         Tz_pi = (Tz_pi - (Nslice_gpu+Nfill0+2)).*s_pi + (Nslice_gpu+Nfill0+2);
        % h_z (skip)
%         toc;
        
% toc;
% tic;
        % interp target (channel)
        t_chn_pos = Eta./delta_d + midchannel_gpu;
        % index and alpha
%         t_chn_floor = floor(t_chn_pos);
%         t_chn_alpha = t_chn_pos - t_chn_floor;
        
              
% toc;
% tic;
        % data iview
        data_iview = gpuArray(dataflow.rawdata(:, :, [iview iview+Nviewprot/2], ishot));
        
        % interp chn 0  
% tic;
%         data_0(:, Nfill0+1:end-Nfill0) = data_iview(t_chn_floor(:), :, 1).*(1-t_chn_alpha(:)) + ...
%             data_iview(t_chn_floor(:)+1, :, 1).*t_chn_alpha(:);
% toc;
% tic;
        data_0(:, Nfill0+1:end-Nfill0) = interp1(channelindex, data_iview(:,:,1), t_chn_pos(:), 'linear', 0);
% toc;
% tic;

%         if ishot == 1
%             data_0(:, 1:Nfill0) = repmat(data_0(: ,Nfill0+1), 1, Nfill0);
%         else
%             data_0(:, 1:Nfill0) = 0;
%         end
%         if ishot == Nshot
%             data_0(:, end-Nfill0+1:end) = repmat(data_0(:, end-Nfill0), 1, Nfill0);
%         else
%             data_0(:, end-Nfill0+1:end) = 0;
%         end

        data_0(:, 1:Nfill0) = repmat(data_0(: ,Nfill0+1), 1, Nfill0).*(ishot == 1);
        data_0(:, end-Nfill0+1:end) = repmat(data_0(:, end-Nfill0), 1, Nfill0).*(ishot == Nshot);
        
% toc;        
% tic;
        % interp on Z
        % self
%         [t_odd, t_even, gamma] = omiga4interp(Tz_self_0, [gamma_coeff1 gamma_coeff2]);

        Tz_floor = floor(Tz_self_0(:));
        s_odd = mod(Tz_floor, 2);
        alpha = Tz_self_0(:) - Tz_floor;
%         beta = sqrt(alpha+1.0);
%         beta = interp1(alpha_intp, t_intp, alpha);
        beta = interp1(index_intp, t_intp, alpha.*Nintp_gpu+1);
        % I know beta = interp1(t_intp, alpha.*Nintp_gpu+1); but which lay in a matlab bug
        t_odd = reshape((Tz_floor+s_odd)./2 + beta(:,2).*s_odd + beta(:,1).*(1-s_odd), imagesize_gpu^2, []);
        t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), imagesize_gpu^2, []);
% toc;
% tic;
        data1 = interp2(data_0(:, 1:2:end), t_odd, index_img, 'linear', 0) + ...
                interp2(data_0(:, 2:2:end), t_even, index_img, 'linear', 0);
        data1 = data1 + interp2(conv2(data_0, [-1 2 -1], 'same'), reshape(Tz_self_0, imagesize_gpu^2, []), index_img, 'linear', 0) ...
                .*reshape(beta(:,3), imagesize_gpu^2, []);
toc;
1;
% % %         t_z_floor = floor(Tz_0(s_0));
% % %         t_z_alpha = Tz_0(s_0) - t_z_floor;
% % %         index_img = reshape(1:imagesize*imagesize, imagesize, imagesize);
%         t_z_floor = floor(Tz_self_0);
%         t_z_alpha = Tz_self_0 - t_z_floor;
% %         
% % %         toc;
% %         % #6
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
% %         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
% %         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
% %         t_z_index = index_img(s_0) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
%         t_z_index = index_img(:) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
% %         t_z_coeff(:,1) = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
%         t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
%                       (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
%                       (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
%                       (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
% %         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
% toc;
% 1;
        % add to image
%         img_shot(s_0) = img_shot(s_0) + sum(data_0(t_z_index).*t_z_coeff, 2);
%         img_shot = img_shot + reshape(sum(data_0(t_z_index).*t_z_coeff, 2), imagesize_gpu, imagesize_gpu, Nextslice);
%         
%         % -1
%         t_z_index = index_img(:) + double(t_z_floor(:)-2).*double(imagesize_gpu)^2;
%         t_z_coeff = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 0
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 1
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 2
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (t_z_alpha(:).*(1-gamma(:))+beta(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
% %         toc;
% 
%         % pi
% %         tic;
%         % interp channel
%         t_chn_alpha = -Eta./delta_d + midchannel_gpu;
%         t_chn_floor = floor(t_chn_alpha);
%         t_chn_alpha = t_chn_alpha - t_chn_floor;
%         
%         % interp chn pi
% %         data_iview(:) = dataflow.rawdata(:, :, iview+Nviewprot/2, ishot);
%         data_0(:, Nfill0+1:end-Nfill0) = data_iview(t_chn_floor(:), :, 2).*(1-t_chn_alpha(:)) + ...
%             data_iview(t_chn_floor(:)+1, :, 2).*t_chn_alpha(:);
%         if ishot == 1
%             data_0(:, 1:Nfill0) = repmat(data_0(: ,Nfill0+1), 1, Nfill0);
%         else
%             data_0(:, 1:Nfill0) = 0;
%         end
%         if ishot == Nshot
%             data_0(:, end-Nfill0+1:end) = repmat(data_0(:, end-Nfill0), 1, Nfill0);
%         else
%             data_0(:, end-Nfill0+1:end) = 0;
%         end
%         
%         
%         t_z_floor = floor(Tz_pi);
%         t_z_alpha = Tz_pi - t_z_floor;
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
% % %         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
% % %         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
% %         t_z_index = index_img(:) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
% % %         t_z_coeff(:,1) = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
% %         t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
% %                       (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
% %                       (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
% %                       (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
% % %         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
% % %         toc;
% % 
% %         % add to image
% % %         img_shot(s_pi) = img_shot(s_pi) + sum(data_0(t_z_index).*t_z_coeff, 2);
% %         img_shot = img_shot + reshape(sum(data_0(t_z_index).*t_z_coeff, 2), imagesize_gpu, imagesize_gpu, Nextslice);
%         
%         % -1
%         t_z_index = index_img(:) + double(t_z_floor(:)-2).*double(imagesize_gpu)^2;
%         t_z_coeff = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 0
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 1
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
%         % 2
%         t_z_index = t_z_index + double(imagesize_gpu)^2;
%         t_z_coeff = (t_z_alpha(:).*(1-gamma(:))+beta(:))./4;
%         img_shot = img_shot + reshape(data_0(t_z_index).*t_z_coeff, imagesize_gpu, imagesize_gpu, Nextslice);
        
%         toc;
    end
    % get img
%     img(:,:,imageindex) = img(:,:,imageindex) + gather(img_shot(:,:,gatherindex));
    fprintf('\n');
toc;
end

img = img.*(pi/Nviewprot/2);

