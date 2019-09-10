function img = backproj2D_1(p, theta, ctrIdx, hond, N, interp)
% back projection 
% img = backproj2D_1(p, theta, ctrIdx, d, N, interp)

% Define the x & y axes for the reconstructed image so that the origin
% (center) is in the spot which RADON would choose.
% center = floor((N + 1)/2);
center = (N + 1)/2;
xleft = -center + 1;
x = (1:N) - 1 + xleft;
x = repmat(x, N, 1);

ytop = center - 1;
y = (N:-1:1).' - N + ytop;
y = repmat(y, 1, N);

interp_method = sprintf('*%s',interp); % Add asterisk to assert
% even-spacing of taxis

% Generate trignometric tables
costheta = cos(theta);
sintheta = sin(theta);

Np = size(p, 1);

switch interp
    case 'linear'
        % Allocate memory for the image
        img = zeros(N,'like',p);
        % interp
        for iview=1:length(theta)
            proj = [p(:,iview); 0; nan];
            t = (x.*costheta(iview) + y.*sintheta(iview)).*hond + ctrIdx;
            t_index = floor(t);
            t_alpha = t - t_index;
            t_index(t_index==0) = Np+1;
            t_index(t_index<0) = Np+2;
            t_index(t_index>Np) = Np+1;
            projContrib = proj(t_index).*(1-t_alpha) + proj(t_index+1).*t_alpha;
            img = img + reshape(projContrib,N,N);
        end
        
    case {'spline','pchip','cubic','v5cubic'}
        % Allocate memory for the image
        img = zeros(N,'like',p);
        % interp
        taxis = ((1:size(p,1)) - ctrIdx)./hond;
        for iview=1:length(theta)
            proj = p(:,iview);
            t = x.*costheta(iview) + y.*sintheta(iview);
            projContrib = interp1(taxis,proj,t(:),interp_method);
            img = img + reshape(projContrib,N,N);
        end
end

img = img.*(pi/length(theta));

return