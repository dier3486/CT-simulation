function img = backproj3Dslope(rawdata, viewangle, XY, Neighb, Chninterp, Zinterp, shotflag, interpmethod, blocksize, GPUonoff)
% back projection 3D axial, GPU free

% I know the size of rawdata is
[Npixel, Nslice, Nviewprot] = size(rawdata);
imagesize2 = size(XY, 1);
% I know imagesize2 = imagesize^2
Nblock = imagesize2/blocksize;


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
    index_img = gpuArray(repmat(single(1:blocksize)', 1, Nextslice));
    Ztarget = repmat(gpuArray(single(z_target)), blocksize, 1);
    data_iview = zeros(Npixel, Nslice, 'single', 'gpuArray');
    data_1 = zeros(blocksize, Nfillslice, 'single', 'gpuArray');
    data_2 = zeros(blocksize, Nextslice, 'single', 'gpuArray');
    blockindex = gpuArray(reshape(single(1:imagesize2), blocksize, Nblock));
else
    img = zeros(imagesize2, Nextslice, 'single');
    index_img = repmat(single(1:blocksize)', 1, Nextslice);
    Ztarget = repmat(single(z_target), blocksize, 1);
    data_iview = zeros(Npixel, Nslice, 'single');
    data_1 = zeros(blocksize, Nfillslice, 'single');
    data_2 = zeros(blocksize, Nextslice, 'single');
    blockindex = reshape(single(1:imagesize2), blocksize, Nblock);
end

costheta = cos(viewangle);
sintheta = sin(viewangle);

for iview = 1:Nviewprot
    % get data per view
    if GPUonoff
        data_iview  = gpuArray(rawdata(:, :, iview));
    else
        data_iview = dataflow.rawdata(:, :, iview);
    end
    for iblock = 1:Nblock
        % sin cos of view angle
        sintheta_iview = sintheta(iview);
        costheta_iview = costheta(iview);
        % X-Y to Zeta-Eta
        Eta = -XY(blockindex(:, iblock),1).*sintheta_iview + XY(blockindex(:, iblock),2).*costheta_iview;
        Zeta = XY(blockindex(:, iblock),1).*costheta_iview + XY(blockindex(:, iblock),2).*sintheta_iview;
        
        % interp on channel direction
        t_chn = Eta./Chninterp.delta_d + Chninterp.midchannel;
        data_1(:, Nfill0+1:end-Nfill0) = interp1(Chninterp.channelindex, data_iview, t_chn(:), 'linear', 0);
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
                t_odd  = reshape((Tz_floor+s_odd)./2 + beta(:,1).*(1-s_odd) + beta(:,2).*s_odd, blocksize, []);
                t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), blocksize, []);
                data_2 = interp2(data_1(:, 1:2:end), t_odd, index_img, 'linear', 0)./2 + ...
                    interp2(data_1(:, 2:2:end), t_even, index_img, 'linear', 0)./2;
                data_2 = data_2 + interp2(conv2(data_1, Zinterp.convL, 'same'), Tz, index_img, 'linear', 0) ...
                    .*reshape(beta(:,3), blocksize, []);
                % I know Zinterpsamp.convL = [-1 2 -1]
            otherwise
                error('Unknown interplation method %s!', interpmethod);
        end
        % add to image
        img(blockindex(:, iblock), :) = img(blockindex(:, iblock), :) + data_2;
    end
end

end