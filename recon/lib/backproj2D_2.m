function img = backproj2D_2(p, theta, ctrIdx, hond, N, interp, centerond)
% back projection 
% img = backproj2D_2(p, theta, ctrIdx, hond, N, interp, centerond)
% p is the rawdata; theta is the viewangle+pi/2; ctrIdx is the midchannel index
% hond is the h/d, where h is the voxel size, d is the ray step length;
% N is the imagesize; interp is interp method;
% centerond is the (centerX, centerY)/d, (centerX, centerY) is the rotation center xy on image coordinates;

% Define the x & y axes for the reconstructed image
[x, y] = ndgrid(-(N-1)/2 : (N-1)/2);
x = x(:).*hond - centerond(:, 1)';
y = y(:).*hond - centerond(:, 2)';

% Generate trignometric tables
costheta = cos(theta);
sintheta = sin(theta);

sizep = size(p);
Np = sizep(1);
if length(sizep) > 2
    Nslice = sizep(2);
else
    Nslice = 1;
    p = reshape(p, sizep(1), 1, sizep(2));
end

if size(x, 2)==1 && Nslice>1
    x = repmat(x, 1, Nslice);
    y = repmat(y, 1, Nslice);
end

% Allocate memory for the image
img = zeros(N,N, Nslice, 'like',p);
switch interp
    case 'linear'
        % interp
        for iview=1:length(theta)
            % fill 0 and nan
            proj = [p(:, :, iview); zeros(1, Nslice); nan(1, Nslice)];
            % projection sample
            t = (-x.*sintheta(iview) + y.*costheta(iview)) + ctrIdx;
            % interpolation index and alpha
            tindex = floor(t);
            t_alpha = t - tindex;
            % index+1
            tindex_p1 = tindex + 1;
            % edge
            tindex(tindex==0) = Np+1;
            tindex(tindex<0) = Np+2;
            tindex(tindex>Np) = Np+1;
            tindex_p1(tindex_p1<=0) = Np+2;
            tindex_p1(tindex_p1>Np+1) = Np+2;
            % slice shift
            tindex = tindex + (Np+2).*(0 : Nslice-1);
            tindex_p1 = tindex_p1 + (Np+2).*(0 : Nslice-1);
            % get sample value by interpolation
            projContrib = proj(tindex(:)).*(1-t_alpha(:)) + proj(tindex_p1(:)).*t_alpha(:);
            % add to image
            img = img + reshape(projContrib, N, N, Nslice);
        end
    otherwise
        1;
        % not support yet
        
%     case {'spline','pchip','cubic','v5cubic'}
%         % interp
%         taxis = ((1:size(p,1)) - ctrIdx)./hond;
%         for iview=1:length(theta)
%             proj = p(:, :, iview);
%             t = x.*costheta(iview) + y.*sintheta(iview);
%             t = t(:);
%             projContrib = interp1(taxis,proj,t(:),interp_method);
%             img = img + reshape(projContrib,N,N);
%         end
end

img = img.*(pi/length(theta)/2);

end