% BP test code
% for 3D Axial no tilt, shots fill up
% half shot

load('E:\data\simulation\TM\test\vol_test1.mat');

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

gpuDevice;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);


% ini image
img = zeros(imagesize, imagesize, Nimage, 'single');

xygrid = gpuArray(single((-(imagesize-1)/2 : (imagesize-1)/2).*h));
[X, Y] = ndgrid(xygrid);

zgrid = gpuArray(single((-(Nslice-1)/2 : (Nslice-1)/2)));
% slicegrid = (-(Nslice-1)/2 : (Nslice-1)/2).*delta_z;
midslice = (Nslice+1)/2;

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% ini GPU buffer
Nslice_gpu = gpuArray(Nslice);
midchannel_gpu = gpuArray(midchannel);
SID_gpu = gpuArray(SID);
imagesize_gpu = gpuArray(imagesize);
delta_z_norm = gpuArray(delta_z/imageincrement/SID);
% Eta = zeros(imagesize, imagesize, 'single', 'gpuArray');
% Zeta = zeros(imagesize, imagesize, 'single', 'gpuArray');
% D = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_self = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_neib = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% gap = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% kg_self = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% kg_next = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% kg_prev = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% s_self = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_prev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_next = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapprev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapnext = false(imagesize, imagesize, Nslice, 'gpuArray');
% t_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% h_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_chn = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv = zeros(imagesize, imagesize, 'single', 'gpuArray');
% data0 = zeros(imagesize, imagesize, Nslice*2, 'single', 'gpuArray');
% data_iview = zeros(Npixel, Nslice*2, 'single', 'gpuArray');
% t_chn_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chn_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_floor = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_chninv_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
% t_z_floor = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_z_alpha = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% beta = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% gamma = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% t_z_index = zeros(imagesize^2*Nslice, 4, 'double', 'gpuArray');
% t_z_coeff = zeros(imagesize^2*Nslice, 4, 'single', 'gpuArray');
index_img = gpuArray(repmat((1:imagesize*imagesize)', Nslice, 1));
img_shot = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');

% coeff
gamma_coeff1 = 0.6;
gamma_coeff2 = 1.5;
gamma_coeff1 = gpuArray(gamma_coeff1);
gamma_coeff2 = gpuArray(gamma_coeff2);


for ishot = 0:Nshot
% for ishot = 1:1
tic;
    if ishot==0
        % first
        imageindex = 1:Nslice/2;
        viewangle = viewangle_prot + startviewangle(1);
        viewshift = 0;
    elseif ishot==Nshot
        % last
        imageindex = (1:Nslice/2) + (Nshot*Nslice-Nslice/2);
        viewangle = viewangle_prot + startviewangle(Nshot);
        viewshift = 0;
    else
        imageindex = (1:Nslice) + (ishot*Nslice-Nslice/2);
        viewangle = viewangle_prot + startviewangle(ishot);
        viewshift = round((startviewangle(ishot) - startviewangle(ishot+1))/(pi*2)*Nviewprot);
    end

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
        t_0 = repmat((D + Zeta).*delta_z_norm, 1, 1, Nslice_gpu);
        t_pi = repmat((D - Zeta).*delta_z_norm, 1, 1, Nslice_gpu);
        gap = Nslice_gpu - (t_self+t_neib).*(Nslice_gpu-1)./2;
        
%         toc;
        % samples from '0' and 'pi' scan
        kg_0 = reshape(zgrid,1,1,[])./t_0 + 1/2;
        kg_pi = reshape(zgrid,1,1,[])./t_pi + 1/2;
        
        % TBC
        
        % to find out then on self or previous or next shot, or in 'gap'
        s_self = (kg_self > -Nslice_gpu/2+1) & (kg_self < Nslice_gpu/2);
        s_prev = (kg_prev < -Nslice_gpu/2) & (kg_self <= -Nslice_gpu/2+1);
        s_next = (kg_next > Nslice_gpu/2+1) & (kg_self >= Nslice_gpu/2);
        s_gapprev = (kg_self <= -Nslice_gpu/2+1) & ~s_prev;
        s_gapnext = (kg_self >= Nslice_gpu/2) & ~s_next;
        
%         toc;
        % interp target (Z)
%         t_z(:) = 0;
%         t_z(s_self) = kg_self(s_self);
%         t_z(s_next) = kg_next(s_next);
%         t_z(s_prev) = kg_prev(s_prev);        
%         t_z(s_gapprev) = -Nslice_gpu/2+1 + (kg_self(s_gapprev) + Nslice_gpu/2 -1).*t_self(s_gapprev)./gap(s_gapprev);
%         t_z(s_gapnext) = Nslice_gpu/2 + (kg_self(s_gapnext) - Nslice_gpu/2).*t_self(s_gapnext)./gap(s_gapnext);        
        
        t_z = kg_self.*s_self + kg_next.*s_next + kg_prev.*s_prev + ...
            (-Nslice_gpu/2+1 + (kg_self + Nslice_gpu/2 -1).*t_self./gap).*s_gapprev + ...
            (Nslice_gpu/2 + (kg_self - Nslice_gpu/2).*t_self./gap).*s_gapnext;
        % + Nslice
        t_z = t_z + Nslice_gpu;
%         toc;
        
        % interp step (Z)
%         h_z(:) = 0;
%         h_z(s_self) = t_self(s_self);
%         h_z(s_prev | s_next) = t_neib(s_prev | s_next);
%         h_z(s_gapprev | s_gapnext) = gap(s_gapprev | s_gapnext);
        
        h_z = t_self.*s_self + t_neib.*(s_prev | s_next) + gap.*(s_gapprev | s_gapnext);
        
%         toc;
        % interp target (channel)
        t_chn = Eta./delta_d + midchannel_gpu;
        t_chninv = -Eta./delta_d + midchannel_gpu;
        % index and alpha
        t_chn_floor = floor(t_chn);
        t_chn_alpha = t_chn - t_chn_floor;
        t_chninv_floor = floor(t_chninv);
        t_chninv_alpha = t_chninv - t_chninv_floor;
              
%         toc;
        % interpolation on channel data
%         data0(:) = 0;
%         % self
%         data0(:, :, Nslice_gpu/2+1:Nslice_gpu*3/2) = reshape(dataflow.rawdata(t_chn_floor(:),:,iview,ishot).*(1-t_chn_alpha(:)) + ...
%             dataflow.rawdata(t_chn_floor(:)+1,:,iview,ishot).*t_chn_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu);
%         toc;
%         % neib
%         if ishot>1
%             iview_prev = mod(viewshift_prev + Nviewprot/2 + iview - 1, Nviewprot) + 1;
%             data0(:, :, 1:Nslice_gpu/2) = ...
%                 reshape(dataflow.rawdata(t_chninv_floor(:), Nslice_gpu/2+1:end, iview_prev, ishot-1).*(1-t_chninv_alpha(:)) + ...
%                 dataflow.rawdata(t_chninv_floor(:)+1, Nslice_gpu/2+1:end, iview_prev, ishot-1).*t_chninv_alpha(:), ...
%                 imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
%         else
%             iview_prev = iview;
%             data0(:, :, 1:Nslice_gpu/2) = repmat(data0(:, :, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
%         end
%         if ishot<Nshot
%             iview_next = mod(viewshift_next + Nviewprot/2 + iview - 1, Nviewprot) + 1;
%             data0(:, :, Nslice_gpu*3/2+1:end) = ...
%                 reshape(dataflow.rawdata(t_chninv_floor(:), 1:Nslice_gpu/2, iview_next, ishot+1).*(1-t_chninv_alpha(:)) + ...
%                 dataflow.rawdata(t_chninv_floor(:)+1, 1:Nslice_gpu/2, iview_next, ishot+1).*t_chninv_alpha(:), ...
%                 imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
%         else
%             iview_next = iview;
%             data0(:, :, Nslice_gpu*3/2+1:end) = repmat(data0(:, :, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
%         end
%         
%         toc;
        
        % data0 #2
%         tic;
        % self
        data_iview(:, Nslice_gpu/2+1:Nslice_gpu*3/2) = dataflow.rawdata(:, :, iview, ishot);
        % neib
        if ishot>1
            iview_prev = mod(viewshift_prev + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data_iview(:, 1:Nslice_gpu/2) = dataflow.rawdata(:, Nslice_gpu/2+1:end, iview_prev, ishot-1);
%         else
%             data_iview(:, 1:Nslice_gpu/2) = repmat(data_iview(:, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
        end
        if ishot<Nshot
            iview_next = mod(viewshift_next + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data_iview(:, Nslice_gpu*3/2+1:end) = dataflow.rawdata(:, 1:Nslice_gpu/2, iview_next, ishot+1);
%         else
%             data_iview(:, Nslice_gpu*3/2+1:end) = repmat(data_iview(:, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
        end
        % interp chn 
        data0(:, :, Nslice_gpu/2+1:Nslice_gpu*3/2) = reshape(data_iview(t_chn_floor(:),Nslice_gpu/2+1:Nslice_gpu*3/2).*(1-t_chn_alpha(:)) + ...
            data_iview(t_chn_floor(:)+1,Nslice_gpu/2+1:Nslice_gpu*3/2).*t_chn_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu);
        % neib
        if ishot>1
            data0(:, :, 1:Nslice_gpu/2) = reshape(data_iview(t_chninv_floor(:), 1:Nslice_gpu/2).*(1-t_chninv_alpha(:)) + ...
                data_iview(t_chninv_floor(:)+1, 1:Nslice_gpu/2).*t_chninv_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
        else
            data0(:, :, 1:Nslice_gpu/2) = repmat(data0(:, :, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
        end
        if ishot<Nshot
            data0(:, :, Nslice_gpu*3/2+1:Nslice_gpu*2) = reshape(data_iview(t_chninv_floor(:), Nslice_gpu*3/2+1:Nslice_gpu*2).*(1-t_chninv_alpha(:)) + ...
                data_iview(t_chninv_floor(:)+1, Nslice_gpu*3/2+1:Nslice_gpu*2).*t_chninv_alpha(:), imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
        else
            data0(:, :, Nslice_gpu*3/2+1:Nslice_gpu*2) = repmat(data0(:, :, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
        end
        
%         toc;
        
        % interp on Z
        t_z_floor = floor(t_z);
        t_z_alpha = t_z - t_z_floor;
%         index_img = reshape(1:imagesize*imagesize, imagesize, imagesize);
        
%         toc;
        % #6
        beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
        gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
        t_z_index = index_img(:) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
%         t_z_coeff(:,1) = ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4;
        t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
                      (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
                      (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
                      (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
%         toc;

        % add to image
        img_shot = img_shot + reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize_gpu, imagesize_gpu, Nslice_gpu);
%         toc;
    end
    % get img
    img(:,:,imageindex) = gather(img_shot);
    fprintf('\n');
toc;
end

img = img.*(pi/Nviewprot/2);

