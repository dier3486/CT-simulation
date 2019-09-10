function D = parallelprojinimage(parallelbeam, Cimage, method)
% parallel projection on image(s)
% D = parallelprojinimage(Np, delta_d, mid_chn, Cimage, viewangle, couch, method)
% shall be faster than the genernal function intersection in as
% intersection.m

if nargin<3
    method = '2D lineinsection';
elseif isempty(method)
    method = '2D lineinsection';
end

% get parameters
Np = parallelbeam.Np;
delta_d = parallelbeam.delta_d;
if isfield(parallelbeam, 'midchannel')
    mid_chn = parallelbeam.midchannel;
else
    mid_chn = (Np+1)/2;
end
if isfield(parallelbeam, 'h')
    h = parallelbeam.h;
else
    h = 1.0;
end
if isfield(parallelbeam, 'viewangle')
    viewangle = parallelbeam.viewangle;
else
    viewangle = 0;
end
% image size
[Nx, Ny, Nz] = size(Cimage);
% viewnumber
Nview = size(viewangle(:), 1);

switch method
    case '2D lineinsection'
        % grid
        Xgrid = (-Nx/2:Nx/2).*h;
        Ygrid = (-Ny/2:Ny/2).*h; 
        % fill a 0 after Cimage
        Cimage = [Cimage(:); 0];
        % d is the distance from ISO to beams
        d = ((1:Np)'-mid_chn).*delta_d;
        % ini D
        D = zeros(Np, Nview, class(Cimage));
        % loop the views 
        for iview = 1:Nview
            theta = repmat(viewangle(iview), Np, 1);
            [dt, Vindex] = linesinimage2D(theta, d, [], 1, Xgrid, Ygrid);
            D(:, iview) = sum(dt.*Cimage(Vindex), 2);
        end
    case '2D linearinterp'
        % fill a 0 after Cimage
        Cimage = [Cimage(:); 0];
        % d_h is the distance from ISO to beams and /h
        d_h = ((1:Np)'-mid_chn).*delta_d./h;
        % ini D
        D = zeros(Np, Nview, class(Cimage));
        % loop the views 
        for iview = 1:Nview
            [inter_a, index_1, index_2, cs_view] = ...
                parallellinearinterp2D(Nx, Ny, d_h, viewangle(iview), Nx*Ny);
            D(:, iview) = sum(Cimage(index_1).*(1-inter_a) + ...
                Cimage(index_2).*inter_a, 2).*(abs(cs_view)*h);
        end
    otherwise
        D = [];
end

return