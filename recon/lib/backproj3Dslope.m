function img = backproj3Dslope(rawdata, viewangle, XY, Neighb, Chninterp, Zinterp, shotflag, interpmethod, viewblock, GPUonoff)
% back projection 3D axial, GPU free

% I know the size of rawdata is
[Npixel, Nslice, Nview] = size(rawdata);
imagesize2 = size(XY, 1);
% I know imagesize2 = imagesize^2
Nviewblock = ceil(Nview/viewblock);

% edge expand 
Nfill0 = Zinterp.Nfill0;
Nfillslice = Nslice + Nfill0*2;

switch shotflag
    case 1
        % first shot
        Nextslice = Nslice + Neighb;
        z_target = -Nslice/2+1 : Nslice/2+Neighb;
    case 2
        % last shot
        Nextslice = Nslice + Neighb;
        z_target = -Nslice/2-Neighb+1 : Nslice/2;
    case 3 
        % only one shot
        Nextslice = Nslice;
        z_target = -Nslice/2+1 : Nslice/2;
    otherwise
        % middle shots
        Nextslice = Nslice + Neighb*2;
        z_target = -Nextslice/2+1 : Nextslice/2;
end

if GPUonoff
    img = zeros(imagesize2, Nextslice, 'single', 'gpuArray');
    index_img = gpuArray(repmat(single(1:imagesize2)', 1, Nextslice));
    Ztarget = repmat(gpuArray(single(z_target)), imagesize2, 1);
    data_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');
    data_1 = zeros(imagesize2, Nfillslice, 'single', 'gpuArray');
    data_2 = zeros(imagesize2, Nextslice, 'single', 'gpuArray');
else
    img = zeros(imagesize2, Nextslice, 'single');
    index_img = repmat(single(1:imagesize2)', 1, Nextslice);
    Ztarget = repmat(single(z_target), imagesize2, 1);
    data_iview = zeros(Npixel, Nslice, 'single');
    data_1 = zeros(imagesize2, Nfillslice, 'single');
    data_2 = zeros(imagesize2, Nextslice, 'single');
end

for iblk = 1:Nviewblock
    if iblk<Nviewblock
        viewindex = (1:viewblock) + (iblk-1)*viewblock;
        Nviewperblk = viewblock;
    else  % iblk == Nviewblock
        viewindex = (iblk-1)*viewblock+1 : Nview;
        Nviewperblk = Nview - viewblock*(Nviewblock-1);
    end
    % get data per block
    if GPUonoff
        datablk = gpuArray(rawdata(:, :, viewindex));
        Nviewperblk = gpuArray(Nviewperblk);
        viewangleblk = gpuArray(viewangle(viewindex));
    else
        datablk = rawdata(:, :, viewindex);
        viewangleblk = viewangle(viewindex);
    end

    % will do
%     % up sampling
%     if upsampling
%         datablk = doubleup(reshape(datablk, Npixel, []), upsampgamma);
%     end
% 
%     % filter
%     if isFilter
%         datablk_f = fft(reshape(datablk, Npixel_up, []), Hlen);
%         datablk_f = ifft(datablk_f.*filter, 'symmetric');
%         datablk = reshape(datablk_f(1:Npixel_up, :), Npixel_up, Nslice, Nviewperblk);
%     end

    costheta = cos(viewangleblk);
    sintheta = sin(viewangleblk);

    % loop the views
    for iview = 1:Nviewperblk
        % X-Y to Zeta-Eta
        Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
        Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);
        % interp on channel direction
        t_chn = Eta./Chninterp.delta_d + Chninterp.midchannel;
        data_1(:, Nfill0+1:end-Nfill0) = interp1(Chninterp.channelindex, datablk(:, :, iview), t_chn(:), 'linear', 0);
        % start and last shot fillup
        switch shotflag
            case 1
                % first shot
                data_1(:, 1:Nfill0) = repmat(data_1(:, Nfill0+1), 1, Nfill0);
            case 2
                % last shot
                data_1(:, end-Nfill0+1:end) = repmat(data_1(:, end-Nfill0), 1, Nfill0);
            case 3
                % only one shot
                data_1(:, 1:Nfill0) = repmat(data_1(: ,Nfill0+1), 1, Nfill0);
                data_1(:, end-Nfill0+1:end) = repmat(data_1(:, end-Nfill0), 1, Nfill0);
            otherwise
                % do nothing
                1;
        end
        % Warning: the branches in a looping could lay in low performance, at least in matlab, though here the switch-case are
        %          determined before the looping to avoid real branches.

        % Zinterp index
        Tz = interp3(Zinterp.Zeta, Zinterp.Eta, Zinterp.zz, Zinterp.t, ...
            repmat(Zeta, 1, Nextslice), repmat(Eta, 1, Nextslice), Ztarget);
        % interp on Z
        switch lower(interpmethod)
            case 'linear'
                data_2 = interp2(data_1, Tz, index_img, 'linear', 0);
            case '4points'
                Tz_floor = floor(Tz(:));
                s_odd = mod(Tz_floor, 2);
                alpha = Tz(:) - Tz_floor;
                % beta & gamma
                beta = interp1(Zinterp.fourpointindex, Zinterp.fourpoint, alpha.*Zinterp.Nfourp+1);
                % 4-points-interp
                t_odd  = reshape((Tz_floor+s_odd)./2 + beta(:,1).*(1-s_odd) + beta(:,2).*s_odd, imagesize2, []);
                t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), imagesize2, []);
                data_2 = interp2(data_1(:, 1:2:end), t_odd, index_img, 'linear', 0)./2 + ...
                    interp2(data_1(:, 2:2:end), t_even, index_img, 'linear', 0)./2;
                data_2 = data_2 + interp2(conv2(data_1, Zinterp.convL, 'same'), Tz, index_img, 'linear', 0) ...
                    .*reshape(beta(:,3), imagesize2, []);
                % I know Zinterpsamp.convL = [-1 2 -1]
            otherwise
                error('Unknown interplation method %s!', interpmethod);
        end
        % add to image
%         img(blockindex(:, iblock), :) = img(blockindex(:, iblock), :) + data_2;
        img = img + data_2;
    end
end

% reshape
img = reshape(img, imagesize2, Nextslice);

end