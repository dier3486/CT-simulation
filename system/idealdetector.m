function detector = idealdetector(detector, direction)
% return ideal detector posisiton
% detector = idealdetector(detector)

if nargin < 1
    detector = [];
end

if nargin < 2
    direction = true;
end

if isempty(detector)
    % system parameters (sample)
    detector.SID = 550;
    detector.SDD = 1000;
    detector.hx_ISO = 0.55;
    detector.hz_ISO = 0.55;
    detector.Npixel = 950;
    detector.Nslice = 16;
    detector.mid_U = 475.25;     % satrt from 1
end

% get detectors position
alpha_1 = asin(detector.hx_ISO/detector.SID);
if direction
    alpha_pixel = ((1:detector.Npixel)-detector.mid_U).*alpha_1 + pi/2;
else
    alpha_pixel = -((1:detector.Npixel)-detector.mid_U).*alpha_1 + pi/2;
end
x_pos = cos(alpha_pixel).*detector.SDD;
y_pos = sin(alpha_pixel).*detector.SDD - detector.SID;
z_pos = ((1:detector.Nslice)-detector.Nslice/2-1/2) .* ...
    (detector.hz_ISO*detector.SDD/detector.SID);

detector.alpha_pixel = alpha_pixel;
detector.position = [repmat(x_pos(:), detector.Nslice, 1), ...
    repmat(y_pos(:), detector.Nslice, 1), repelem(z_pos(:), detector.Npixel, 1)];

return

