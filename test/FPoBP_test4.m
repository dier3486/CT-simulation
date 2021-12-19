% test data
% load F:\data-Dier.Z\3.2Head338\test\data_bnh_1122.mat

status.nodename='Boneharden';

% prm
imagesize = single(prmflow.recon.imagesize);
h = single(prmflow.recon.FOV/prmflow.recon.imagesize);
delta_d = single(prmflow.recon.delta_d);
hond = h/delta_d;
Npixel = prmflow.recon.Npixel;
midchannel = prmflow.recon.midchannel;
Nviewprot = prmflow.recon.Nviewprot;
reconcenter = prmflow.recon.center;
delta_view = prmflow.recon.delta_view;
startviewangle = prmflow.recon.startviewangle;
% Nimage = prmflow.recon.Nimage;
Nslice = single(prmflow.recon.Nslice);
FOV = prmflow.recon.FOV;
effFOV = single(min(FOV*1.2, prmflow.recon.maxFOV));
effNp = min(floor(effFOV/delta_d) + 1, Npixel);
SID = prmflow.recon.SID;
Neighb = single(prmflow.recon.Neighb);
Nextslice = prmflow.recon.Nextslice;
Nedge = single(2);
% Nshot = prmflow.recon.Nshot;

Nshot = 3;
Nimage = Nshot*Nslice;

subview = 4;

% bone
% calibration table
bonecorr = prmflow.corrtable.(status.nodename);
bonecurve = reshape(bonecorr.bonecurve, bonecorr.bonecurvelength, []);

% enhance bone edge
if isfield(prmflow.pipe.(status.nodename), 'edgeenhance')
    edgeenhance = prmflow.pipe.(status.nodename).edgeenhance;
else
    edgeenhance = true;
end
if isfield(prmflow.pipe.(status.nodename), 'edgekernel')
    edgekernel = prmflow.pipe.(status.nodename).edgekernel;
else
    edgekernel = 2.0;
end
if isfield(prmflow.pipe.(status.nodename), 'edgescale')
    edgescale = prmflow.pipe.(status.nodename).edgescale;
else
    edgescale = 1.0;
end
if edgeenhance
    minvalue = min(bonecurve(:,1)) - 20;
    pixelsize = prmflow.recon.FOV/prmflow.recon.imagesize;
    ImageBE = bonecorrEdgeEnhance(dataflow.image(:,:,1:Nimage), minvalue, pixelsize, edgekernel, edgescale);
else
    ImageBE = dataflow.image(:,:,1:Nimage);
end
% use bone curve to adjust image
BoneImage = GetBoneImg(ImageBE, bonecurve);

image_out = zeros(imagesize, imagesize, Nimage, 'single');

% prepare
% view angle
viewangle = mod((0:Nviewprot/2-1).*delta_view + startviewangle(1) + pi/2, pi*2);
viewangle = reshape(viewangle, subview, []);
Nview_perit = size(viewangle, 2);

eta_C = reconcenter(1).*sin(viewangle) - reconcenter(2).*cos(viewangle);
indexstart_p = floor(midchannel + (-effFOV/2 + eta_C)./delta_d);
indexstart_p(indexstart_p<1) = 1;
indexstart_n = ceil(midchannel*2 - indexstart_p);
indexstart_n(indexstart_n>Npixel) = Npixel;

viewangle = gpuArray(viewangle);

% ctrIdx = midchannel+eta_C./delta_d+1-indexstart;
channelpos = ((1:Npixel)'-midchannel).*delta_d;
channelindex = gpuArray(single(1:effNp)');
maxR = effFOV/2/delta_d;
channelpos_h = gpuArray(channelpos./h);
reconcenter_h = gpuArray(reconcenter./h);

% XY
xygrid = single((-(imagesize-1)/2 : (imagesize-1)/2).*h);
[X, Y] = ndgrid(xygrid);
Sxy = gpuArray(sqrt(X.^2+Y.^2) <= maxR*delta_d);
XY = gpuArray([X(Sxy) Y(Sxy)] - reconcenter);
%  XY = gpuArray([X(:) Y(:)] - reconcenter);
Nxy = size(XY, 1);
SID_h = gpuArray(SID/h);
% zz_samp = single(-Nslice+1:Nslice);
Zgrid_bp = gpuArray((1:Nslice)-1/2);

% tic;
% BBH table
HCscale = 1000;
Dscale = gpuArray(1/HCscale/bonecorr.curvescale(1));
bonescale = gpuArray(1/HCscale/bonecorr.curvescale(2));
bonecorr.order = reshape(bonecorr.order, 1, []);
curvematrix = gpuArray(reshape(bonecorr.curvematrix, bonecorr.order));
efffilter = interp1( bonecorr.beamposition, bonecorr.effbeamfilter, channelpos, 'linear', 'extrap');
efffilter = gpuArray(efffilter./bonecorr.curvescale(3));
mubonmuw = gpuArray(single(bonecorr.refrencebonemu/bonecorr.refrencemu));
Nw = 400; Nb = 200; Nf = 20;
gridW = gpuArray(single(linspace(0, 0.5, Nw)));
gridB = gpuArray(single(linspace(0, 2.0, Nb)));
gridF = gpuArray(single(linspace(0.1, 1.0, Nf)));
[gBB, gWW, gFF] = meshgrid(gridB, gridW, gridF);
BBHmatrix = polyval3dm(curvematrix, gWW, gBB, gFF).*mubonmuw-1;
% toc;
% Filter
filter = gpuArray(prmflow.recon.filter);
Hlen = length(filter);

% costheta = cos(viewangle);
% sintheta = sin(viewangle);

% BV
tau_BV = 3.0;

mu_BV0 = 0.33;
mu_adap = 0.15;
lambda = 0.01;
mu_F = 1.0;
mu_L = 0.03/mu_F;
% filters
Klr = gpuArray(filterdesign('ram-lak', Npixel, delta_d, 2.0));
% main filter
alpha_lr = 0;
alpha_filt = 2.0;
Kfilt = gpuArray(filterdesign('hann', Npixel, delta_d, alpha_filt));
Kfilt = Kfilt.*(1-alpha_lr) + Klr.*alpha_lr;
% noise filter
alpha_filt = 2.0;
alpha_lr = 0.9;
Kfilt_BV = filterdesign('hann', Npixel, delta_d, alpha_filt);
% Kfilt_BV = gpuArray(filterdesign('ram-lak', pbm.Npixel, pbm.delta_d, alpha_filt));
Kfilt_BV = Kfilt_BV.*(1-alpha_lr) + Klr.*alpha_lr;

% F^2
Kf2 = Kfilt.*Kfilt_BV.*pi;
% or
Kf1 = Kfilt_BV.*pi;

% iteration ini
image0 = gpuArray(dataflow.image(:,:,1:Nimage)).*Sxy;
% mu ini
Gu = BregmanTV2D(image0, mu_BV0, lambda, []);
mu_BV = mu_BV0./(abs(image0-Gu).*mu_adap+1);

% iteration 
Niter = 1;
alpha_iter = 0.25;
image1 = image0;
Rerr = zeros(Niter*subview, Nimage);
% 1st step BV
Gu = BregmanTV2D(image1, mu_BV, lambda, Gu);
for iiter = 1 : Niter*subview
    % new mu_BV
    mu_BV = mu_BV0./(abs(image1-Gu).*mu_adap+1);
    % left vector 
    u = (image1+Gu)./2 + 1.0i.*(image1-Gu).*mu_F;
    
    isub = mod(iiter-1, subview) + 1;
    % loop shots
    image_diff = reshape(image0.*0, imagesize^2, Nimage);
    for ishot = 0 : Nshot
        switch ishot
            case 0
                % FP
                img_index = [ones(1, Neighb) 1:Nslice/2+Nedge];
                Nslice_ishot = Nslice/2;
                Zgrid_fp = gpuArray((-Nslice/2+1:0) - 1/2);
                Zshift_fp = ones(1, Nslice/2, 'single', 'gpuArray').*(Nslice/2) + Neighb + 1/2;
                % BP
                Zgrid_bp = gpuArray((Nslice/2+1:Nslice)-1/2);
                Q1 = zeros(Nxy, Nslice/2+Nedge+Neighb, 'single', 'gpuArray');
                Q2 = zeros(Nxy, Nslice/2+Nedge+Neighb, 'single', 'gpuArray');
                indexQ1 = gpuArray([ones(1, Neighb)  (1:Nslice/2)  ones(1, Nedge).*(Nslice/2)]);
                indexQ2 = gpuArray([ones(1, Neighb)+Nslice/2  (1:Nslice/2)+Nslice/2  ones(1, Nedge).*Nslice]);
                interpZ_shift = Neighb-Nslice/2+1/2;
                imgbk_index = 1:Nslice/2;
            case Nshot
                % FP
                img_index = [(-Nslice/2-Nedge+1 : 0)+Nslice*Nshot  ones(1, Neighb).*(Nshot*Nslice)];
                Nslice_ishot = Nslice/2;
                Zgrid_fp = gpuArray((1:Nslice/2) - 1/2);
                Zshift_fp = zeros(1, Nslice/2, 'single', 'gpuArray') + Nedge + 1/2;
                % BP
                Zgrid_bp = gpuArray((1:Nslice/2)-1/2);
                Q1 = zeros(Nxy, Nslice/2+Nedge+Neighb, 'single', 'gpuArray');
                Q2 = zeros(Nxy, Nslice/2+Nedge+Neighb, 'single', 'gpuArray');
                indexQ1 = gpuArray([ones(1, Nedge)  (1:Nslice/2)  ones(1, Neighb).*(Nslice/2)]);
                indexQ2 = gpuArray([ones(1, Nedge)+Nslice/2  (1:Nslice/2)+Nslice/2  ones(1, Neighb).*Nslice]);
                interpZ_shift = Nedge+1/2;
                imgbk_index = (-Nslice/2+1 : 0) + Nslice*Nshot;
            otherwise
                % FP
                img_index = (-Nslice/2-Nedge+1 : Nslice/2+Nedge) + Nslice*ishot;
                Nslice_ishot = Nslice;
                Zgrid_fp = gpuArray([1:Nslice/2 -Nslice/2+1:0] - 1/2);
                Zshift_fp = [zeros(1, Nslice/2, 'single', 'gpuArray') ones(1, Nslice/2, 'single', 'gpuArray').*Nslice] + Nedge + 1/2;
                % BP
                Zgrid_bp = gpuArray((1:Nslice)-1/2);
                Q1 = zeros(Nxy, Nslice+Nedge*2, 'single', 'gpuArray');
                Q2 = zeros(Nxy, Nslice+Nedge*2, 'single', 'gpuArray');
                indexQ1 = gpuArray([ones(1, Nedge)  (1:Nslice/2)  (1:Nslice/2)+Nslice*3/2  ones(1, Nedge).*(2*Nslice)]);
                indexQ2 = gpuArray([ones(1, Nedge)+Nslice  (1:Nslice/2)+Nslice  (1:Nslice/2)+Nslice/2  ones(1, Nedge).*Nslice]);
                interpZ_shift = Nedge+1/2;
                imgbk_index = (-Nslice/2+1 : Nslice/2) + Nslice*ishot;
        end

        if iiter == 1
            imagebone_ishot = gpuArray(BoneImage(:,:,img_index)).*Sxy;
        end
        
        index_voxel = gpuArray(repmat(single(1:Nxy)', 1, Nslice_ishot));
        image_fix = zeros(Nxy, Nslice_ishot, 'single', 'gpuArray');
        % image_fix = zeros(imagesize, imagesize, Nslice, 'single', 'gpuArray');
        
        tic;
        for iview = 1:Nview_perit
            %     tic;
            theta = viewangle(isub, iview);
            % sintheta_iview = sintheta(iview);
            % costheta_iview = costheta(iview);
            
            % FP
            dh_iview = [channelpos_h(indexstart_p(isub,iview):indexstart_p(isub,iview)+effNp-1); ...
                       -channelpos_h(indexstart_n(isub,iview)-effNp+1:indexstart_n(isub,iview))];
            % interpXY
            [interpX, interpY, cs_view] = parallellinearinterp2D2(imagesize, imagesize, dh_iview, theta, reconcenter_h);
            interpY_rep = repmat(interpY, 1, 1, Nslice_ishot);
            interpX_rep = repmat(interpX, 1, 1, Nslice_ishot);
            % move XY to center
            interpX = interpX - (imagesize+1)/2 - reconcenter_h(1);
            interpY = interpY - (imagesize+1)/2 - reconcenter_h(2);
            % interpZ
            Eta = -interpX.*sin(theta) + interpY.*cos(theta);
            Zeta = interpX.*cos(theta) + interpY.*sin(theta);
            Zeta(effNp+1:end, :) = -Zeta(effNp+1:end, :);
            Zscale = sqrt(SID_h.^2 - Eta.^2);
            Zscale = (Zscale + Zeta)./Zscale;
            interpZ = Zscale(:)*Zgrid_fp + Zshift_fp;
            interpZ = reshape(interpZ, effNp*2, imagesize, Nslice_ishot);
            % project
            P0 = sum(interp3(u(:,:,img_index), interpY_rep, interpX_rep, interpZ, 'linear', 0), 2).*(abs(cs_view)*h);
            P0 = reshape(P0, effNp, Nslice_ishot*2);
            % F^2
            Dv = ifft(fft(real(P0), Hlen).*Kfilt + fft(imag(P0), Hlen).*Kf2, 'symmetric');
            % bone
            if iiter == 1
                DB = sum(interp3(imagebone_ishot, interpY_rep, interpX_rep, interpZ, 'linear', 0), 2).*(abs(cs_view)*h);
                DB = reshape(DB, effNp, Nslice_ishot*2);
                % remove negative
                D0 = real(P0);
                D0 = D0.*(D0>0);
                DB = DB.*(DB>0);
                
                effF_ivew = repmat([efffilter(indexstart_p(isub,iview):indexstart_p(isub,iview)+effNp-1) ...
                                    efffilter(indexstart_n(isub,iview)-effNp+1:indexstart_n(isub,iview))], 1, Nslice_ishot);
                Dfix = interp3(gBB, gWW, gFF, BBHmatrix, DB.*bonescale, D0.*Dscale, effF_ivew).*DB;
                Df = ifft(fft(Dfix, Hlen).*filter, 'symmetric');
                D = Df(1:effNp, :) + 1.0i.*Dv(1:effNp, :);
            else
                D = Dv(1:effNp, :);
            end

            % BP
            % X-Y to Zeta-Eta
            Eta = -XY(:,1).*sin(theta) + XY(:,2).*cos(theta);
            Zeta = XY(:,1).*cos(theta) + XY(:,2).*sin(theta);
            
            t_chn1 = Eta./delta_d + midchannel + 1 - indexstart_p(isub,iview);
            %         P1 = interp1(channelindex, Df(:, 1:2:end), t_chn1(:), 'linear', 0);
            t_chn2 = -Eta./delta_d + midchannel + effNp - indexstart_n(isub,iview);
            %         P2 = interp1(channelindex, Df(:, 2:2:end), t_chn2(:), 'linear', 0);
            P12 = [interp1(channelindex, D(:, 1:2:end), t_chn1(:), 'linear', 0)  interp1(channelindex, D(:, 2:2:end), t_chn2(:), 'linear', 0)];
            
            SetaD = sqrt(SID.^2 - Eta.^2);
            
            Zscale_p = SetaD./(SetaD+Zeta);
            Zscale_n = SetaD./(SetaD-Zeta);
            Z1_p = Zscale_p(:)*Zgrid_bp;
            Z1_n = Zscale_n(:)*Zgrid_bp;
            s1_p = Z1_p<=Nslice/2;
            s1_n = Z1_n<=Nslice/2;
            Z2_p = Zscale_p(:)*(Zgrid_bp-Nslice) + Nslice;
            Z2_n = Zscale_n(:)*(Zgrid_bp-Nslice) + Nslice;
            
            Q1 = P12(:, indexQ1);
            Q2 = P12(:, indexQ2);
            
            interpZ1 = Z1_p.*s1_p + Z2_n.*~s1_p + interpZ_shift;
            interpZ2 = Z1_n.*s1_n + Z2_p.*~s1_n + interpZ_shift;
            
            %     image_fix = image_fix + reshape(interp2(Q1, interpZ1, index_img, 'linear'), imagesize, imagesize, Nslice) + ...;
            %              reshape(interp2(Q2, interpZ2, index_img, 'linear'), imagesize, imagesize, Nslice);
            image_fix = image_fix + interp2(Q1, interpZ1, index_voxel, 'linear') + ...;
                interp2(Q2, interpZ2, index_voxel, 'linear');
            %     toc;
        end
        image_fix = image_fix.*(pi/Nview_perit/4);
        image_diff(Sxy, imgbk_index) = image_fix;
        toc;
    end
    image_diff = reshape(image_diff, imagesize, imagesize, Nimage);
    
    % add to origine image
    if iiter == 1
        %1 BH
        image1 = image0 + real(image_diff);
        %2 noise
        r = image0-imag(image_diff);
        
        % update image0
        image0 = image0 + real(image_diff);
        Gu = Gu + real(image_diff);
    else
        r = image0-image_diff;
    end
    Rerr(iiter, :) = gather(sqrt(sum(reshape(r, [], Nimage).^2, 1))./imagesize^2);
    image1 = image1 + r.*alpha_iter;
    % BV
    Gu = BregmanTV2D(image1, mu_BV, lambda, Gu);
    % I know the result is (image1+Gu)./2
end



function ImgOut = GetBoneImg(ImgIn, BoneCurve)
minValue = min(BoneCurve(:,1));
maxValue = max(BoneCurve(:,1));
ImgIn(ImgIn < minValue) = 0;
ImgIn(ImgIn > maxValue) = maxValue;
idx = find(ImgIn > 0);
ImgIn(idx)=interp1(BoneCurve(:,1), BoneCurve(:,2), ImgIn(idx));
ImgOut = ImgIn;
end