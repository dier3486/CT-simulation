function [interpX, interpY, interpZ, cs_view] = parallellinearinterp3D3(imagesize, d_h, Za, Zb, viewangle, centerXY_h)
% the subfunction to calculation the sample points' position in 3D parallel forward projection.
% [interpX, interpY, interpZ, cs_view] = ...
%   parallellinearinterp3D3(imagesize, d_h, Za, Zb, viewangle, centerXY_h);
% then the projection reads
% P = sum(interp3(image, interpX, interpY, interpZ, 'linear', 0), 2).*abs(cs_view);
% where
%   for homocentric, Za = Zslice*sqrt(SID_h^2-d_h^2)/SID_h - Phi*Cd,  Zb = Zslice/SID_h;
%   for concyclic, Za = Zslice - Phi*Cd,  Zb = Zslice/sqrt(SID_h^2-d_h^2);
% 

if nargin<6
    centerXY_h = [0 0];
end

gridmax = -(max(imagesize)-1)/2 : (max(imagesize)-1)/2;
C = centerXY_h*[-sin(viewangle); cos(viewangle)];
Np = size(d_h, 1);
Nz = size(Za, 2);
Nrepb = Np/size(Zb, 1);

tantheta = tan(viewangle);
phase_flag = abs(tantheta) > 1;

% phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);
if phase_flag
    % y-like
    cs_view = -csc(viewangle);
    cottheta = 1/tantheta;

    interpX = (gridmax.*cottheta + (imagesize(1)+1)/2) + repmat((d_h + C).*cs_view, Nz, 1);
    interpY = repmat(gridmax + (imagesize(2)+1)/2, Np*Nz, 1);
    interpZ = repelem((-Zb(:).*cs_view)*gridmax, Nrepb, 1) + reshape(-Zb.*(d_h.*cottheta - centerXY_h(2)*cs_view) + Za, [], 1);

else
    % x-like
    cs_view = sec(viewangle);

    interpX = repmat(gridmax + (imagesize(1)+1)/2, Np*Nz, 1);
    interpY = (gridmax.*tantheta + (imagesize(2)+1)/2) + repmat((d_h + C).*cs_view, Nz, 1);
    interpZ = repelem((Zb(:).*cs_view)*gridmax, Nrepb, 1) + reshape(Zb.*(d_h.*tantheta - centerXY_h(1)*cs_view) + Za, [], 1);

end

end