% BP test code
% for 3D Axial no tilt, shots fill up

% inputs are dataflow, prmflow
if exist('df0', 'var')
    dataflow = df0;
    clear df0;
end

if exist('pf0', 'var')
    prmflow = pf0;
    clear pf0;
end

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
% imageincrement = prmflow.recon.imageincrement;
imageincrement = delta_z;
imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
% set center to (0, 0)
imagecenter(:, 1:2) = 0;
% no tilt

h = FOV/imagesize;


% ini image
img = zeros(imagesize, imagesize, Nimage);

xygrid = (-(imagesize-1)/2 : (imagesize-1)/2).*h;
[X, Y] = ndgrid(xygrid);

zgrid = (-(Nslice-1)/2 : (Nslice-1)/2);
% slicegrid = (-(Nslice-1)/2 : (Nslice-1)/2).*delta_z;
midslice = (Nslice+1)/2;

for ishot = 1:Nshot
% for ishot = 1:1
    imageindex = (1:Nslice) + (ishot-1)*Nslice;
    viewangle = startviewangle(ishot) + linspace(0, pi*2, Nviewprot+1);
    viewangle = viewangle(1:end-1);
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
    for iview = 1:Nviewprot
%         tic;
        % .
        fprintf('.');
        % X-Y to Eta-Zeta
        Eta = -X.*sintheta(iview) + Y.*costheta(iview);
        Zeta = X.*costheta(iview) + Y.*sintheta(iview);
        
        % samples step on Z
        D = sqrt(SID^2 - Eta.^2);
        t_self = (D + Zeta).*(delta_z/imageincrement/SID);
        t_neib = (D - Zeta).*(delta_z/imageincrement/SID);
        gap = Nslice - (t_self+t_neib).*(Nslice-1)./2;
        
        % samples from 'self' and next/previous shots
        kg0 = reshape(zgrid,1,1,[])./t_self + 1/2;
        kg_next = (reshape(zgrid,1,1,[]) - Nslice)./t_neib + 1/2 + Nslice;
        kg_prev = (reshape(zgrid,1,1,[]) + Nslice)./t_neib + 1/2 - Nslice;
        % to find out then on self or previous or next shot, or in 'gap'
        s1 = kg0 > -Nslice/2+1;
        s2 = kg0 < Nslice/2;
        s0 = s1 & s2;
        s_prev = (kg_prev < -Nslice/2) & ~s1;
        s_next = (kg_next > Nslice/2+1) & ~s2;
        s_gap1 = ~s1 & ~s_prev;
        s_gap2 = ~s2 & ~s_next;
               
        % interp target (Z)
        t_z = zeros(size(kg0));
        t_z(s0) = kg0(s0);
        t_z(s_next) = kg_next(s_next);
        t_z(s_prev) = kg_prev(s_prev);
        t_self = repmat(t_self, 1, 1, Nslice);
        t_neib = repmat(t_neib, 1, 1, Nslice);
        gap = repmat(gap, 1, 1, Nslice);
        t_z(s_gap1) = -Nslice/2+1 + (kg0(s_gap1) + Nslice/2 -1).*t_self(s_gap1)./gap(s_gap1);
        t_z(s_gap2) = Nslice/2 + (kg0(s_gap2) - Nslice/2).*t_self(s_gap2)./gap(s_gap2);
        % + Nslice
        t_z = t_z + Nslice;
        
        % interp step (Z)
        h_z = zeros(size(kg0));
        h_z(s0) = t_self(s0);
        h_z(s_prev | s_next) = t_neib(s_prev | s_next);
        h_z(s_gap1 | s_gap2) = gap(s_gap1 | s_gap2);
        
        % interp target (channel)
        t_chn = Eta./delta_d + midchannel;
        t_chninv = -Eta./delta_d + midchannel;
        
        % interpolation on channel data
        data0 = zeros(imagesize, imagesize, Nslice*2);
        % self
        t_chn_floor = floor(t_chn);
        t_chn_alpha = t_chn - t_chn_floor;
        data0(:, :, Nslice/2+1:Nslice*3/2) = reshape(dataflow.rawdata(t_chn_floor(:),:,iview,ishot).*(1-t_chn_alpha(:)) + ...
            dataflow.rawdata(t_chn_floor(:)+1,:,iview,ishot).*t_chn_alpha(:), imagesize, imagesize, Nslice);
        % neib
        t_chninv_floor = floor(t_chninv);
        t_chninv_alpha = t_chninv - t_chninv_floor;
        if ishot>1
            iview_prev = mod(viewshift_prev + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data0(:, :, 1:Nslice/2) = ...
                reshape(dataflow.rawdata(t_chninv_floor(:), Nslice/2+1:end, iview_prev, ishot-1).*(1-t_chninv_alpha(:)) + ...
                dataflow.rawdata(t_chninv_floor(:)+1, Nslice/2+1:end, iview_prev, ishot-1).*t_chninv_alpha(:), ...
                imagesize, imagesize, Nslice/2);
        else
            iview_prev = iview;
            data0(:, :, 1:Nslice/2) = repmat(data0(:, :, Nslice/2+1), 1, 1, Nslice/2);
        end
        if ishot<Nshot
            iview_next = mod(viewshift_next + Nviewprot/2 + iview - 1, Nviewprot) + 1;
            data0(:, :, Nslice*3/2+1:end) = ...
                reshape(dataflow.rawdata(t_chninv_floor(:), 1:Nslice/2, iview_next, ishot+1).*(1-t_chninv_alpha(:)) + ...
                dataflow.rawdata(t_chninv_floor(:)+1, 1:Nslice/2, iview_next, ishot+1).*t_chninv_alpha(:), ...
                imagesize, imagesize, Nslice/2);
        else
            iview_next = iview;
            data0(:, :, Nslice*3/2+1:end) = repmat(data0(:, :, Nslice*3/2), 1, 1, Nslice/2);
        end
        
        % interp on Z
        t_z_floor = floor(t_z);
        t_z_alpha = t_z - t_z_floor;
        index_img = reshape(1:imagesize*imagesize, imagesize, imagesize);
        
%         % #1
%         t_z_index1 = repmat(index_img, 1, 1, Nslice) + (t_z_floor-1).*imagesize^2;
%         t_z_index2 = repmat(index_img, 1, 1, Nslice) + t_z_floor.*imagesize^2;
%         data1 = data0(t_z_index1).*(1-t_z_alpha) + data0(t_z_index2).*t_z_alpha;

%         % #2
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
%         tmp1 = (1 - t_z_alpha(:) + beta(:)).*h2on1(:);
%         tmp2 = (t_z_alpha(:) + beta(:)).*h2on3(:); 
%         t_z_coeff = [tmp1  3-t_z_alpha(:).*2-tmp1  1+t_z_alpha(:).*2-tmp2  tmp2]./4;
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);

%         % #3
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
%         gamma_max = 1.0;
%         gamma_scale = 0.5;
%         gamma = h_z.*gamma_scale;
% %         gamma = (1./h_z).*gamma_scale;
% %         gamma = 1.0 - abs(h_z-1.0);
%         gamma(gamma<0) = 0;
%         gamma(gamma>gamma_max) = gamma_max;
%         tmp1 = ((1-t_z_alpha(:)+beta(:))./4 + (-2+t_z_alpha(:)+beta(:)).*gamma(:)./8).*h2on1(:);
%         tmp2 = ((t_z_alpha(:)+beta(:))./4 + (-1-t_z_alpha(:)+beta(:)).*gamma(:)./8).*h2on3(:);
%         t_z_coeff = [tmp1  (3-t_z_alpha(:).*2)./4+(1-t_z_alpha(:).*2).*gamma(:)./8-tmp1 ...
%                      (1+t_z_alpha(:).*2)./4+(-1+t_z_alpha(:).*2).*gamma(:)./8-tmp2  tmp2];
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
        
%         % #4
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
%         gamma_max = 1.6;
%         gamma_scale = 0.8;
%         gamma = h_z.*gamma_scale;
% %         gamma = (1./h_z).*gamma_scale;
% %         gamma = 1.0 - abs(h_z-1.0);
%         gamma(gamma<0) = 0;
%         gamma(gamma>gamma_max) = gamma_max;
%         tmp1 = ((1-t_z_alpha(:)+beta(:))./4 + (-1+t_z_alpha(:)).*gamma(:)./4).*h2on1(:);
%         tmp2 = ((t_z_alpha(:)+beta(:))./4 + (-t_z_alpha(:)).*gamma(:)./4).*h2on3(:);
%         t_z_coeff = [tmp1  (3-t_z_alpha(:).*2)./4+(1-t_z_alpha(:).*2).*gamma(:)./4-tmp1 ...
%                      (1+t_z_alpha(:).*2)./4+(-1+t_z_alpha(:).*2).*gamma(:)./4-tmp2  tmp2];
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
        
%         % #5
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
% %         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
% %         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
%         gamma_max = 1.0;
%         gamma_scale = 0.5;
%         gamma = h_z.*gamma_scale;
% % %         gamma = (1./h_z).*gamma_scale;
% % %         gamma = 1.0 - abs(h_z-1.0);
% %         gamma(gamma<0) = 0;
%         gamma(gamma>gamma_max) = gamma_max;
%         t_z_coeff = [ (1-t_z_alpha(:)+beta(:))./4+(-1+t_z_alpha(:)-beta(:)).*gamma(:)./4 ...
%                       (2-t_z_alpha(:)-beta(:))./4+(2-t_z_alpha(:).*3+beta(:)).*gamma(:)./4 ...
%                       (1+t_z_alpha(:)-beta(:))./4+(-1+t_z_alpha(:).*3+beta(:)).*gamma(:)./4 ...
%                       (t_z_alpha(:)+beta(:))./4+(-t_z_alpha(:)-beta(:)).*gamma(:)./4];        
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
        
        % #6
        beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
        gamma_c1 = 0.6;
        gamma_c2 = 1.5;
        gamma = gamma_c1./sqrt(1-t_z_alpha.*(1-t_z_alpha).*gamma_c2);
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
        t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
        t_z_coeff = [ ((1-t_z_alpha(:)).*(1-gamma(:))+beta(:))./4 ...
                      (2-t_z_alpha(:)-beta(:)+(2-t_z_alpha(:).*3).*gamma(:))./4 ...
                      (1+t_z_alpha(:)-beta(:)+(t_z_alpha(:).*3-1).*gamma(:))./4 ...
                      (t_z_alpha(:).*(1-gamma(:))+beta(:))./4];     
        data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);
        
%         % #7
%         beta = 1/2 - sqrt(1/4+t_z_alpha.*(1-t_z_alpha));
%         h2on1 = cat(3, ones(imagesize, imagesize), h_z(:, :, 2:end)./h_z(:, :, 1:end-1));
%         h2on3 = cat(3, h_z(:, :, 1:end-1)./h_z(:, :, 2:end), ones(imagesize, imagesize));
%         t_z_index = repmat(index_img(:), Nslice, 1) + (t_z_floor(:) + [-2 -1 0 1]).*imagesize^2;
%         gamma_scale = 1.0;
%         gamma = gamma_scale./sqrt(1-t_z_alpha.*(1-t_z_alpha));
%         
%         tmp1 = ((1-t_z_alpha(:)+beta(:))./4 + (-1+t_z_alpha(:)).*gamma(:)./4).*h2on1(:);
%         tmp2 = ((t_z_alpha(:)+beta(:))./4 + (-t_z_alpha(:)).*gamma(:)./4).*h2on3(:);
%         t_z_coeff = [tmp1  (3-t_z_alpha(:).*2)./4+(1-t_z_alpha(:).*2).*gamma(:)./4-tmp1 ...
%             (1+t_z_alpha(:).*2)./4+(-1+t_z_alpha(:).*2).*gamma(:)./4-tmp2  tmp2];
%         data1 = reshape(sum(data0(t_z_index).*t_z_coeff, 2), imagesize, imagesize, Nslice);

        % add to image
        img(:, :, imageindex) = img(:, :, imageindex) + data1;
%         toc;
    end
    fprintf('\n');
end

img = img.*(pi/Nviewprot/2);

