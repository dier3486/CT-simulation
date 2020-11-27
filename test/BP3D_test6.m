% BP test code
% for 3D Axial no tilt, shots fill up
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
imageincrement = delta_z;
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

% zgrid = gpuArray(single((-(Nslice-1)/2 : (Nslice-1)/2)));
% zgrid = gpuArray(single(1/2 : Nslice-1/2));
Zgrid = gpuArray(single(-Nextslice/2+1/2 : Nextslice/2-1/2));
% slicegrid = (-(Nslice-1)/2 : (Nslice-1)/2).*delta_z;
% midslice = (Nslice+1)/2;
% index_np = [Nextslice/2+1:Nextslice  1:Nextslice/2];
index_np = gpuArray([(Nslice+1:Nextslice)  repmat(Nextslice, 1, Nslice-Nextslice/2)  ...
            ones(1,Nslice-Nextslice/2)  1:Nextslice-Nslice]);

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% edge expand 
Nfill0 = gpuArray(single(4));
Nextslice = gpuArray(Nextslice);

% ini GPU buffer
Nslice_gpu = gpuArray(Nslice);
midchannel_gpu = gpuArray(midchannel);
SID_gpu = gpuArray(SID);
imagesize_gpu = gpuArray(imagesize);
delta_z_norm = gpuArray(delta_z/imageincrement/SID);
Eta = zeros(imagesize, imagesize, 'single', 'gpuArray');
Zeta = zeros(imagesize, imagesize, 'single', 'gpuArray');
D = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_0 = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_pi = zeros(imagesize, imagesize, 'single', 'gpuArray');
gap = zeros(imagesize, imagesize, 'single', 'gpuArray');
kg_0 = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
kg_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% kg_prev = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
s_0 = false(imagesize, imagesize, Nextslice, 'gpuArray');
s_pi = false(imagesize, imagesize, Nextslice, 'gpuArray');
s_gap = false(imagesize, imagesize, Nextslice, 'gpuArray');
% s_gapprev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapnext = false(imagesize, imagesize, Nslice, 'gpuArray');
Tz_0 = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
Tz_pi = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');
% h_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_chn = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv = zeros(imagesize, imagesize, 'single', 'gpuArray');
data_0 = zeros(imagesize, imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
% data_0 = zeros(imagesize, imagesize, Nslice+Nfill0*2, 'single', 'gpuArray');
data_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');
t_chn_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_chn_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_z_floor = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_z_alpha = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% beta = zeros(imagesize, imagesize, Nslice*2, 'single', 'gpuArray');
% gamma = zeros(imagesize, imagesize, Nslice*2, 'single', 'gpuArray');
% t_z_index = zeros(imagesize^2*Nslice, 4, 'double', 'gpuArray');
% t_z_coeff = zeros(imagesize^2*Nslice, 4, 'single', 'gpuArray');
index_img = gpuArray(repmat((1:imagesize*imagesize)', Nextslice, 1));
img_shot = zeros(imagesize, imagesize, Nextslice, 'single', 'gpuArray');

Nshift1 = reshape(repelem([Nslice_gpu/2-1 -Nslice_gpu/2], 1, Nextslice/2), 1, 1, []);
Nshift2 = reshape(repelem([-Nslice_gpu Nslice_gpu], 1, Nextslice/2), 1, 1, []); 
Nshift3 = reshape(repelem([2 Nslice_gpu+Nfill0+2], 1, Nextslice/2), 1, 1, []); 

% coeff
gamma_coeff1 = 0.6;
gamma_coeff2 = 1.5;
gamma_coeff1 = gpuArray(gamma_coeff1);
gamma_coeff2 = gpuArray(gamma_coeff2);

 

for ishot = 1:Nshot
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
        % .
        fprintf('.');
%         tic;
        % X-Y to Eta-Zeta
        Eta = -X.*sintheta(iview) + Y.*costheta(iview);
        Zeta = X.*costheta(iview) + Y.*sintheta(iview);
        % I know the +pi is -Eta and -Zeta
        
        % samples step on Z
        D = sqrt(SID_gpu^2 - Eta.^2);
        t_0 = (D + Zeta).*delta_z_norm;
        t_pi = (D - Zeta).*delta_z_norm;
        gap = Nslice_gpu - (t_0+t_pi).*(Nslice_gpu-1)./2;
        
%         toc;
        % samples from '0' and 'pi' scan
        kg_0 = reshape(Zgrid,1,1,[])./t_0 + 1/2;
        kg_pi = reshape(Zgrid,1,1,[])./t_pi + 1/2;
        % I know kg_neib_0 = -flip(kg_self_0)
        
        s_0 = (kg_0 < Nslice_gpu/2) & (kg_0 > -Nslice_gpu/2+1);   
        s_pi = (kg_pi < Nslice_gpu/2) & (kg_pi > -Nslice_gpu/2+1);
        % I know s_neib_0 = s_self_0(index_np)
        s_gap = ~s_0 & ~s_pi(:,:,index_np);
        % I know s_gap_pi = s_gap_0(index_np)
        
        % t_z
%         Tz_0(:,:,1:Nslice_gpu) = -Nslice_gpu/2;     Tz_0(:,:,Nslice_gpu+1:end) = Nslice_gpu/2;
%         Tz_pi(:,:,1:Nslice_gpu) = -Nslice_gpu/2;     Tz_pi(:,:,Nslice_gpu+1:end) = Nslice_gpu/2;
        Tz_0 = (kg_pi(:,:,index_np) + Nshift2).*s_pi(:,:,index_np) + kg_0.*s_0 + ((kg_0 + Nshift1).*(t_0./gap) - Nshift1).*s_gap + Nslice_gpu/2+Nfill0;
        Tz_pi = (kg_0(:,:,index_np) + Nshift2).*s_0(:,:,index_np) + kg_pi.*s_pi + ((kg_pi + Nshift1).*(t_pi./gap) - Nshift1).*s_gap(:,:,index_np) + Nslice_gpu/2+Nfill0;
%         toc;
        1;
        
%         Tz_0(Tz_0<2) = 2;    Tz_0(Tz_0>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
%         Tz_pi(Tz_pi<2) = 2;    Tz_pi(Tz_pi>Nslice_gpu+Nfill0+2) = Nslice_gpu+Nfill0+2;
        
        s_0 = (Tz_0 >= 2) & (Tz_0 <= (Nslice_gpu+Nfill0+2));
        Tz_0 = Tz_0.*s_0 + Nshift3.*~s_0;
        s_pi = (Tz_pi >= 2) & (Tz_pi <= (Nslice_gpu+Nfill0+2));
        Tz_pi = Tz_pi.*s_pi + Nshift3.*~s_pi;
        
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
        
        
        % interp step (Z)
%         h_z(:) = 0;
%         h_z(s_self) = t_self(s_self);
%         h_z(s_prev | s_next) = t_neib(s_prev | s_next);
%         h_z(s_gapprev | s_gapnext) = gap(s_gapprev | s_gapnext);
        
%         h_z = t_self.*s_self + t_neib.*(s_prev | s_next) + gap.*(s_gapprev | s_gapnext);
        
        tic;
        % interp target (channel)
        t_chn_alpha = Eta./delta_d + midchannel_gpu;
        % index and alpha
        t_chn_floor = floor(t_chn_alpha);
        t_chn_alpha = t_chn_alpha - t_chn_floor;
        
              
%         toc;

        
        % data0 #2
        % interp chn 0
        data_iview(:) = dataflow.rawdata(:, :, iview, ishot);
        data_0(:) = 0;
        data_0(:, :, Nfill0+1:end-Nfill0) = reshape(data_iview(t_chn_floor(:), :).*(1-t_chn_alpha(:)) + ...
            data_iview(t_chn_floor(:)+1, :).*t_chn_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu);
        
%         toc;
        % interp on Z
%         t_z_floor = floor(Tz_0(s_0));
%         t_z_alpha = Tz_0(s_0) - t_z_floor;
%         index_img = reshape(1:imagesize*imagesize, imagesize, imagesize);
        t_z_floor = floor(Tz_0);
        t_z_alpha = Tz_0 - t_z_floor;
        
%         toc;
        % #6
        beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
        gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = index_img(s_0) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
        t_z_index = index_img(:) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
%         t_z_coeff(:,1) = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
        t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
                      (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
                      (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
                      (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
%         toc;

        % add to image
%         img_shot(s_0) = img_shot(s_0) + sum(data_0(t_z_index).*t_z_coeff, 2);
        img_shot = img_shot + reshape(sum(data_0(t_z_index).*t_z_coeff, 2), imagesize_gpu, imagesize_gpu, Nextslice);        
        toc;

        % pi
        tic;
        % interp channel
        t_chn_alpha = -Eta./delta_d + midchannel_gpu;
        t_chn_floor = floor(t_chn_alpha);
        t_chn_alpha = t_chn_alpha - t_chn_floor;
        
        % interp chn pi
        data_iview(:) = dataflow.rawdata(:, :, iview+Nviewprot/2, ishot);
        data_0(:) = 0;
        data_0(:, :, Nfill0+1:end-Nfill0) = reshape(data_iview(t_chn_floor(:), :).*(1-t_chn_alpha(:)) + ...
            data_iview(t_chn_floor(:)+1, :).*t_chn_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu);
        
        t_z_floor = floor(Tz_pi(s_pi));
        t_z_alpha = Tz_pi(s_pi) - t_z_floor;
        beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
        gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
        t_z_index = index_img(s_pi) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
%         t_z_coeff(:,1) = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
        t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
                      (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
                      (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
                      (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
%         toc;

        % add to image
        img_shot(s_pi) = img_shot(s_pi) + sum(data_0(t_z_index).*t_z_coeff, 2);
        toc;
    end
    % get img
    img(:,:,imageindex) = img(:,:,imageindex) + gather(img_shot(:,:,gatherindex));
    fprintf('\n');
toc;
end

img = img.*(pi/Nviewprot/2);

