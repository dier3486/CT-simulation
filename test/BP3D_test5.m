% BP test code
% for 3D Axial tilt, shots fill up

load('E:\data\simulation\TM\test\tiltbp_test1.mat');

% % inputs are dataflow, prmflow
% if exist('df0', 'var')
%     dataflow = df0;
%     clear df0;
% end
% 
% if exist('pf0', 'var')
%     prmflow = pf0;
%     clear pf0;
% end

detector = prmflow.system.detector;
SID = detector.SID;
SDD = detector.SDD;
delta_z = detector.hz_ISO;
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
imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
gantrytilt = single(prmflow.recon.gantrytilt);
imageincrement = prmflow.recon.imageincrement;
% imageincrement = delta_z;
shotcouchstep = prmflow.recon.shotcouchstep;

% gantrytilt = 0;
% imagecenter(:,[1 2]) = 0;
imagecenter0 = mean(imagecenter, 1);
FOV = 300;
h = FOV/imagesize;
imageincrement = imageincrement*cos(gantrytilt);

% I know: -zgrid(:).*delta_z.*sin(gantrytilt) == imagecenter(1:Nslice,2), for thin image film

gpuDevice;

% ini image
img = zeros(imagesize, imagesize, Nimage, 'single');

xygrid = gpuArray(single((-(imagesize-1)/2 : (imagesize-1)/2).*h));
[X0, Y0] = ndgrid(xygrid);
% X, Y, (Z) grid for images
X = X0 + reshape(imagecenter(1:Nslice, 1), 1, 1, Nslice);
Y = Y0 + reshape(imagecenter(1:Nslice, 2), 1, 1, Nslice);
% X, Y, (Z) grid for slices
X0 = X0 + imagecenter0(1);
Y0 = Y0 + imagecenter0(2);

zgrid = gpuArray(single((-(Nslice-1)/2 : (Nslice-1)/2)));
zgrid = reshape(zgrid,1,1,[]);
slicegrid = gpuArray((-(Nslice-1)/2 : (Nslice-1)/2).*(delta_z/SID));
midslice = gpuArray((Nslice+1)/2);
centerstep = gpuArray((-shotcouchstep)*sin(gantrytilt));

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% ini GPU buffer
Nslice_gpu = gpuArray(Nslice);
midchannel_gpu = gpuArray(midchannel);
SID_gpu = gpuArray(SID);
imagesize_gpu = gpuArray(imagesize);
delta_z_norm = gpuArray(delta_z/imageincrement/SID);
Eta = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
Zeta = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
D = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_self = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_neib = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
gap = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
kg_self = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% kg_next = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% kg_prev = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
kg_neib = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
s_self = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_prev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_next = false(imagesize, imagesize, Nslice, 'gpuArray');
s_neib = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapprev = false(imagesize, imagesize, Nslice, 'gpuArray');
% s_gapnext = false(imagesize, imagesize, Nslice, 'gpuArray');
s_gap = false(imagesize, imagesize, Nslice, 'gpuArray');
t_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
% h_z = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_chn = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_chninv = zeros(imagesize, imagesize, 'single', 'gpuArray');

t_chn_index = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_chn_alpha = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_chninv_index = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_chninv_alpha = zeros(imagesize, imagesize, 'single', 'gpuArray');
t_z_index1 = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_z_index2 = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_z_alpha = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
beta = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
gamma = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
t_z_index = zeros(imagesize^2*Nslice, 4, 'double', 'gpuArray');
t_z_coeff = zeros(imagesize^2*Nslice, 4, 'single', 'gpuArray');
index_img = gpuArray(repmat((1:imagesize*imagesize)', Nslice, 1));
img_shot = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
data0 = zeros(imagesize, imagesize, Nslice*2, 'single', 'gpuArray');
data_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');
% dataneib_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');

a_stilt = zeros(1, Nslice, 'single', 'gpuArray');
b1_stilt = zeros(imagesize^2, 1, 'single', 'gpuArray');
b2_stilt = zeros(imagesize^2, Nslice, 'single', 'gpuArray');
c_stilt = (X0(:).^2 + Y0(:).^2 - SID^2);
% w_stilt = zeros(imagesize^2, Nslice, 'single', 'gpuArray');
zNsshift = zeros(1, 1, Nslice, 'single', 'gpuArray');
singantrytilt = gpuArray(sin(gantrytilt));
cosgantrytilt = gpuArray(cos(gantrytilt));

% coeff
gamma_coeff1 = 0.6;
gamma_coeff2 = 1.5;
gamma_coeff1 = gpuArray(gamma_coeff1);
gamma_coeff2 = gpuArray(gamma_coeff2);
Nshot = 1;

for ishot = 1:Nshot
% for ishot = 1:1
tic;
    imageindex = (1:Nslice) + (ishot-1)*Nslice;
    viewangle = viewangle_prot + startviewangle(ishot);
    if ishot>1
        viewshift_prev = round((startviewangle(ishot) - startviewangle(ishot-1))/(pi*2)*Nviewprot);
    else
        viewshift_prev = 0;
    end
    if ishot<Nshot
        viewshift_next = round((startviewangle(ishot) - startviewangle(ishot+1))/(pi*2)*Nviewprot);
    else
        viewshift_next = 0;
    end
%     viewangle = linspace(0, pi*2, Nviewprot);
    costheta = cos(viewangle);
    sintheta = sin(viewangle);
%     Xis = X(:) - imagecenter(sliceindex, 1)';
%     Yis = Y(:) - imagecenter(sliceindex, 2)';
    
    % ini
    img_shot(:) = 0;
    for iview = 1:Nviewprot
        % .
        fprintf('.');
%         tic;
        % X-Y to Eta-Zeta
        Eta = -X.*sintheta(iview) + Y.*costheta(iview);
        Zeta = X.*costheta(iview) + Y.*sintheta(iview);
        
        % samples step on Z
        D = sqrt(SID_gpu^2 - Eta.^2);
        t_self = (D + Zeta).*delta_z_norm;
        % neighbers
        % next
%         Eta_next = Eta + centerstep*costheta(iview);
%         D = centerstep*sintheta(iview) + sqrt(SID_gpu^2 - (Eta + centerstep*costheta(iview)).^2);
        D(:, :, Nslice_gpu/2+1:end) = centerstep*sintheta(iview) + ...
            sqrt(SID_gpu^2 - (Eta(:, :, Nslice_gpu/2+1:end) + centerstep*costheta(iview)).^2);
%         t_next = (D - Zeta).*delta_z_norm;
        % prev
%         Eta_prev = Eta - centerstep*costheta(iview);
%         D = -centerstep*sintheta(iview) + sqrt(SID_gpu^2 - (Eta - centerstep*costheta(iview)).^2);
        D(:, :, 1:Nslice_gpu/2) = -centerstep*sintheta(iview) + ...
            sqrt(SID_gpu^2 - (Eta(:, :, 1:Nslice_gpu/2) - centerstep*costheta(iview)).^2);
%         t_prev = (D - Zeta).*delta_z_norm;
        t_neib = (D - Zeta).*delta_z_norm;
%         toc;
        % edge
%         edge_next = floor((Nslice_gpu - (Nslice_gpu-1).*t_self./2)./t_next + 1/2);
%         edge_next(edge_next>Nslice_gpu/2) = Nslice_gpu/2;
%         edge_prev = floor((Nslice_gpu - (Nslice_gpu-1).*t_self./2)./t_prev + 1/2);
%         edge_prev(edge_prev>Nslice_gpu/2) = Nslice_gpu/2;
        edge = floor((Nslice_gpu - (Nslice_gpu-1).*t_self./2)./t_neib + 1/2);
        edge = edge.*(edge<=Nslice_gpu/2) + (edge>Nslice_gpu/2).*(Nslice_gpu/2);
        % edge(edge>Nslice_gpu/2) = Nslice_gpu/2;
%         t_neib = (D - Zeta).*delta_z_norm;
%         gap_next = Nslice_gpu - t_self.*(Nslice_gpu-1)./2 - t_next.*(edge_next-1/2);
%         gap_prev = Nslice_gpu - t_self.*(Nslice_gpu-1)./2 - t_prev.*(edge_prev-1/2);
        gap = Nslice_gpu - t_self.*(Nslice_gpu-1)./2 - t_neib.*(edge-1/2);
        gap = gap + (gap<eps).*eps;
%         gap(abs(gap)<eps) = eps;
        
        % samples from 'self' and next/previous shots
        kg_self = zgrid./t_self + 1/2;
        kg_neib(:,:,1:Nslice_gpu/2) = (zgrid(:,:,1:Nslice_gpu/2) + Nslice_gpu)./t_neib(:,:,1:Nslice_gpu/2) + 1/2 - Nslice_gpu;
        kg_neib(:,:,Nslice_gpu/2+1:end) = (zgrid(:,:,Nslice_gpu/2+1:end) - Nslice_gpu)./t_neib(:,:,Nslice_gpu/2+1:end) + 1/2 + Nslice_gpu;
%         kg_next = (zgrid - Nslice_gpu)./t_next + 1/2 + Nslice_gpu;
%         kg_prev = (zgrid + Nslice_gpu)./t_prev + 1/2 - Nslice_gpu;
        
        % to find out then on self or previous or next shot, or in 'gap'
        s_self = (kg_self > -Nslice_gpu/2+1) & (kg_self < Nslice_gpu/2);
        s_neib = abs(zgrid) > Nslice_gpu - (edge-1/2).*t_neib;
%         s_prev = zgrid < -Nslice_gpu + (edge_prev-1/2).*t_prev;
%         s_next = zgrid > Nslice_gpu - (edge_next-1/2).*t_next;
        s_gap = ~s_self & ~s_neib;
%         s_gapprev = (kg_self <= -Nslice_gpu/2+1) & ~s_prev;
%         s_gapnext = (kg_self >= Nslice_gpu/2) & ~s_next;
        
%         toc;
        % interp target (Z)
%         t_z(:) = 0;
%         t_z(s_self) = kg_self(s_self);
%         t_z(s_next) = kg_next(s_next);
%         t_z(s_prev) = kg_prev(s_prev);
        
%         t_z = kg_self.*s_self + kg_next.*s_next + kg_prev.*s_prev;
%         t_z(s_gapprev) = -Nslice_gpu/2+1 + (kg_self(s_gapprev) + Nslice_gpu/2 -1).*t_self(s_gapprev)./gap_prev(s_gapprev);
%         t_z(s_gapnext) = Nslice_gpu/2 + (kg_self(s_gapnext) - Nslice_gpu/2).*t_self(s_gapnext)./gap_next(s_gapnext);        
        
%         t_z = kg_self.*s_self + kg_next.*s_next + kg_prev.*s_prev + ...
%             (-Nslice_gpu/2+1 + (kg_self + Nslice_gpu/2 -1).*t_self./gap_prev).*s_gapprev + ...
%             (Nslice_gpu/2 + (kg_self - Nslice_gpu/2).*t_self./gap_next).*s_gapnext;
        t_self = t_self./gap;
        zNsshift(1:Nslice_gpu/2) = -Nslice_gpu/2+1;  zNsshift(Nslice_gpu/2+1:end) = Nslice_gpu/2;
        t_z = kg_self.*s_self + kg_neib.*s_neib + (kg_self.*t_self + zNsshift.*(1-t_self)).*s_gap;
        
%         toc;
        % interp index and alpha (Z)
        t_z_index1 = floor(t_z);
        t_z_alpha = t_z - t_z_index1;
        % + Nslice
        t_z_index1 = t_z_index1 + Nslice_gpu;
        t_z_index2 = t_z_index1 + 1;
        % gap
%         zNsshift(1:Nslice_gpu/2) = 0;  zNsshift(Nslice_gpu/2+1:end) = Nslice_gpu*2+1;
%         edge(:,:,Nslice_gpu/2+1:end) = -edge(:,:,Nslice_gpu/2+1:end);
        t_z_index1(:,:,1:Nslice_gpu/2) = t_z_index1(:,:,1:Nslice_gpu/2).*~s_gap(:,:,1:Nslice_gpu/2) + ...
            edge(:,:,1:Nslice_gpu/2).*s_gap(:,:,1:Nslice_gpu/2);
        t_z_index2(:,:,Nslice_gpu/2+1:end) = t_z_index2(:,:,Nslice_gpu/2+1:end).*~s_gap(:,:,Nslice_gpu/2+1:end) + ...
            (Nslice_gpu*2+1-edge(:,:,Nslice_gpu/2+1:end)).*s_gap(:,:,Nslice_gpu/2+1:end);
%         t_z_index1(s_gapprev) = edge_prev(s_gapprev);
%         t_z_index2(s_gapnext) = Nslice_gpu*2+1 - edge_next(s_gapnext);
        
        % interp step (Z)
%         h_z(:) = 0;
%         h_z(s_self) = t_self(s_self);
%         h_z(s_prev | s_next) = t_neib(s_prev | s_next);
%         h_z(s_gapprev | s_gapnext) = gap(s_gapprev | s_gapnext);
        
%         h_z = t_self.*s_self + t_neib.*(s_prev | s_next) + gap.*(s_gapprev | s_gapnext);
%         h_z = t_self.*s_self + t_prev.*s_prev + t_next.*s_next + gap_prev.*s_gapprev + gap_next.*s_gapnext;
        
%         toc;
        % interp target (channel)
        % tilt shift for slices
        a_stilt = cosgantrytilt^2 + slicegrid.^2.*singantrytilt.^2 - 2.*slicegrid.*singantrytilt.*cosgantrytilt.*sintheta(iview);
        b1_stilt = -(X0(:).*costheta(iview) + Y0(:).*sintheta(iview)).*cosgantrytilt;
        b2_stilt = Y0(:)*(slicegrid.*singantrytilt);
%         b_stilt = b_stilt./((SID*cosgantrytilt)^2 + slicegrid.^2.*singantrytilt.^2 - ...
%             2.*slicegrid.*singantrytilt.*(SID*cosgantrytilt).*sintheta(iview));
        % c_stilt done
        Y = reshape((sqrt((b1_stilt + b2_stilt).^2-c_stilt.*a_stilt)-(b1_stilt + b2_stilt))./a_stilt.*slicegrid, ...
            imagesize, imagesize, Nslice);
%         Y = reshape(Y0(:) + b_stilt.*singantrytilt, imagesize, imagesize, Nslice);
        Y = Y0 + Y.*singantrytilt;
        Eta = -X0.*sintheta(iview) + Y.*costheta(iview);
        t_chn_alpha = Eta./delta_d + midchannel_gpu;

        % tilt shift for slices +pi
        a_stilt = cosgantrytilt^2 + slicegrid.^2.*singantrytilt.^2 + 2.*slicegrid.*singantrytilt.*cosgantrytilt.*sintheta(iview);
%         b_stilt = (X0(:).*costheta(iview) + Y0(:).*sintheta(iview)).*cosgantrytilt + Y0(:)*(slicegrid.*singantrytilt);
%         b_stilt = (sqrt(b_stilt.^2-c_stilt.*a_stilt)-b_stilt)./a_stilt.*slicegrid;
        Y = reshape((sqrt((-b1_stilt + b2_stilt).^2-c_stilt.*a_stilt)-(-b1_stilt + b2_stilt))./a_stilt.*slicegrid, ...
            imagesize, imagesize, Nslice);
        Y = Y0 + Y.*singantrytilt;
        Eta = -X0.*sintheta(iview) + Y.*costheta(iview);
        t_chninv_alpha = -Eta./delta_d + midchannel_gpu;
        
%         toc;
%         fprintf('#2\n');
        
        % index and alpha
        t_chn_index = floor(t_chn_alpha);  t_chn_alpha = t_chn_alpha - t_chn_index; 
        t_chninv_index = floor(t_chninv_alpha);  t_chninv_alpha = t_chninv_alpha - t_chninv_index;
%         toc;
        % outof range?
        t_chn_index(t_chn_index<1) = 1;
        t_chn_index(t_chn_index>Npixel-1) = Npixel-1;
        t_chninv_index(t_chninv_index<1) = 1;
        t_chninv_index(t_chninv_index>Npixel-1) = Npixel-1; 
        % slice shift
        t_chn_index = t_chn_index + reshape((0:Nslice_gpu-1).*Npixel, 1, 1, []);
        t_chninv_index = t_chninv_index + reshape((0:Nslice_gpu-1).*Npixel, 1, 1, []);
              
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
        % self data
        data_iview = gpuArray(dataflow.rawdata(:, :, iview, ishot));
        % interp chn self
        data0(:, :, Nslice_gpu/2+1:Nslice_gpu*3/2) = reshape(data_iview(t_chn_index).*(1-t_chn_alpha) + ...
            data_iview(t_chn_index+1).*t_chn_alpha, imagesize_gpu, imagesize_gpu, Nslice_gpu);
        % neib data
%         data_iview = cat(2, repmat(data_iview(:, end), 1, Nslice_gpu/2), repmat(data_iview(:, 1), 1, Nslice_gpu/2));
        if ishot>1
            iview_prev = mod(viewshift_prev + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data_iview(:, Nslice_gpu/2+1:end) = dataflow.rawdata(:, Nslice_gpu/2+1:end, iview_prev, ishot-1);
%         else
%             data_iview(:, 1:Nslice_gpu/2) = repmat(data_iview(:, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
        end
        if ishot<Nshot
            iview_next = mod(viewshift_next + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data_iview(:, 1:Nslice_gpu/2) = dataflow.rawdata(:, 1:Nslice_gpu/2, iview_next, ishot+1);
%         else
%             data_iview(:, Nslice_gpu*3/2+1:end) = repmat(data_iview(:, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
        end
        % interp chn neib
        data0(:, :, [Nslice_gpu*3/2+1:Nslice_gpu*2 1:Nslice_gpu/2]) = reshape(data_iview(t_chninv_index).*(1-t_chninv_alpha) + ...
                data_iview(t_chninv_index+1).*t_chninv_alpha, imagesize_gpu, imagesize_gpu, Nslice_gpu);
        if ishot == 1
            data0(:, :, 1:Nslice_gpu/2) = repmat(data0(:, :, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
        end
        if ishot==Nshot
            data0(:, :, Nslice_gpu*3/2+1:Nslice_gpu*2) = repmat(data0(:, :, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
        end
        
%         if ishot>1
%             data0(:, :, 1:Nslice_gpu/2) = reshape(data_iview(t_chninv_index(:, :, Nslice_gpu/2+1:end)).*(1-t_chninv_alpha(:, Nslice_gpu/2+1:end)) + ...
%                 data_iview(t_chninv_index(:, Nslice_gpu/2+1:end)+1).*t_chninv_alpha(:, Nslice_gpu/2+1:end), imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
%         else
%             data0(:, :, 1:Nslice_gpu/2) = repmat(data0(:, :, Nslice_gpu/2+1), 1, 1, Nslice_gpu/2);
%         end
%         if ishot<Nshot
%             data0(:, :, Nslice_gpu*3/2+1:Nslice_gpu*2) = reshape(data_iview(t_chninv_index(:, 1:Nslice_gpu/2)).*(1-t_chninv_alpha(:, 1:Nslice_gpu/2)) + ...
%                 data_iview(t_chninv_index(:, 1:Nslice_gpu/2)+1).*t_chninv_alpha(:, 1:Nslice_gpu/2), imagesize_gpu, imagesize_gpu, Nslice_gpu/2);
%         else
%             data0(:, :, Nslice_gpu*3/2+1:Nslice_gpu*2) = repmat(data0(:, :, Nslice_gpu*3/2), 1, 1, Nslice_gpu/2);
%         end
        
%         toc;
        
        % interp on Z

%         index_img = reshape(1:imagesize*imagesize, imagesize, imagesize);
        
%         toc;
        % #6
        beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
        gamma = gamma_coeff1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_coeff2);
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = index_img(:) + (double(t_z_floor(:)) + [-2 -1 0 1]).*double(imagesize_gpu)^2;
        t_z_index = index_img(:) + double([t_z_index1(:) + [-2 -1] t_z_index2(:) + [-1 0]]).*double(imagesize_gpu)^2;
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

