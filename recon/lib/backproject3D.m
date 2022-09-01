function img = backproject3D(img, rawdata, viewangle, XY, Chninterp, Zinterp, shotflag, interpmethod, ...
    index_img, Ztarget, databuff1, databuff2)
% back projection 3D axial, subfunction

% parameters
Nviewperblk = length(viewangle);
imagesize2 = size(XY, 1);
Nfill0 = Zinterp.Nfill0;
Nextslice = size(Ztarget, 2);

% cos and sin of viewangle
costheta = cos(viewangle);
sintheta = sin(viewangle);

% loop the views
for iview = 1:Nviewperblk
    % X-Y to Zeta-Eta
    Eta = -XY(:, 1).*sintheta(iview) + XY(:, 2).*costheta(iview);
    Zeta = XY(:, 1).*costheta(iview) + XY(:, 2).*sintheta(iview);
    % interp on channel direction
    t_chn = Eta./Chninterp.delta_d + Chninterp.midchannel;
    databuff1(:, Nfill0+1:end-Nfill0) = interp1(Chninterp.channelindex, rawdata(:, :, iview), t_chn(:), 'linear', 0);
    % start and last shot fillup
    switch shotflag
        case 1
            % first shot
            databuff1(:, 1:Nfill0) = repmat(databuff1(:, Nfill0+1), 1, Nfill0);
        case 2
            % last shot
            databuff1(:, end-Nfill0+1:end) = repmat(databuff1(:, end-Nfill0), 1, Nfill0);
        case 3
            % only one shot
            databuff1(:, 1:Nfill0) = repmat(databuff1(: ,Nfill0+1), 1, Nfill0);
            databuff1(:, end-Nfill0+1:end) = repmat(databuff1(:, end-Nfill0), 1, Nfill0);
        otherwise
            % do nothing
            1;
    end

    % Zinterp index
    Tz = interp3(Zinterp.Zeta, Zinterp.Eta, Zinterp.zz, Zinterp.t, ...
        repmat(Zeta, 1, Nextslice), repmat(Eta, 1, Nextslice), Ztarget);
    % interp on Z
    switch lower(interpmethod)
        case 'linear'
            databuff2 = interp2(databuff1, Tz, index_img, 'linear', 0);
        case '4points'
            Tz_floor = floor(Tz(:));
            s_odd = mod(Tz_floor, 2);
            alpha = Tz(:) - Tz_floor;
            % beta & gamma
            beta = interp1(Zinterp.fourpointindex, Zinterp.fourpoint, alpha.*Zinterp.Nfourp+1);
            % 4-points-interp
            t_odd  = reshape((Tz_floor+s_odd)./2 + beta(:,1).*(1-s_odd) + beta(:,2).*s_odd, imagesize2, []);
            t_even = reshape((Tz_floor-s_odd)./2 + beta(:,1).*s_odd + beta(:,2).*(1-s_odd), imagesize2, []);
            databuff2 = interp2(databuff1(:, 1:2:end), t_odd, index_img, 'linear', 0)./2 + ...
                interp2(databuff1(:, 2:2:end), t_even, index_img, 'linear', 0)./2;
            databuff2 = databuff2 + interp2(conv2(databuff1, Zinterp.convL, 'same'), Tz, index_img, 'linear', 0) ...
                .*reshape(beta(:,3), imagesize2, []);
            % I know Zinterpsamp.convL = [-1 2 -1]
        otherwise
            error('Unknown interplation method %s!', interpmethod);
    end
    % add to image
    img = img + databuff2;
end



end