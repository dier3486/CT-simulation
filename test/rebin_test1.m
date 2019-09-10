addpath(genpath('D:/matlab/CTsimulation/'));

% inputs
% detector = load('D:/matlab/CTsimulation/system/detectorframe/detectorpos_ideal_1000.mat');

focalspot = [0, -detector.SID, 0];
Npixel = detector.Npixel;

Nview = 1440;
% Nview = 3;
% viewangle = single(linspace(0, pi*2, Nview+1));
viewangle = linspace(0, pi*2, Nview+1);
viewangle = viewangle(1:end-1);

% raw0 = ones(Npixel, Nview);
% raw0 = rand(Npixel, Nview);
% raw0 = reshape(1:Npixel*Nview, Npixel, Nview);
raw0 = D(1:Npixel, :);

% go
% fan angles
y = detector.position(1:Npixel, 2) - focalspot(2);
x = detector.position(1:Npixel, 1) - focalspot(1);
fanangles = atan2(y, x);
% d is the distance from ray to ISO
Lxy = sqrt(x.^2+y.^2);
d1 = (detector.position(1:Npixel, 1).*focalspot(2) - detector.position(1:Npixel, 2).*focalspot(1))./Lxy;
% or
if righttoleft
    d = -detector.SID.*cos(fanangles);
else
    % or negative order
    d = detector.SID.*cos(fanangles);
end

% rebin 1
delta_view = pi*2/Nview;
f = fanangles./delta_view;
viewindex = floor(f);
interalpha = f-viewindex;
viewindex = viewindex + 1;  % start from 0
startvindex = mod(max(viewindex), Nview)+1;
viewindex = repmat(viewindex, 1, Nview) + repmat(0:Nview-1, Npixel, 1);
vindex1 = mod(viewindex-1, Nview).*Npixel + repmat((1:Npixel)', 1, Nview);
vindex2 = mod(viewindex, Nview).*Npixel + repmat((1:Npixel)', 1, Nview);

A = zeros(Npixel, Nview);
A(vindex1) = raw0.*repmat(1-interalpha, 1, Nview);
A(vindex2) = A(vindex2) + raw0.*repmat(interalpha, 1, Nview);

% start angle for first rebin view
A = [A(:, startvindex:end) A(:, 1:startvindex-1)];
viewangle = [viewangle(startvindex:end) viewangle(1:startvindex-1)];
startviewangle = viewangle(1);

% rebin 2
delta_t = detector.hx_ISO;
% I know the d is in negative order
mid_t = mod(detector.mid_U, 1)*delta_t;
t1 = ceil((min(d)+mid_t)/delta_t);
t2 = floor((max(d)+mid_t)/delta_t);
Nreb = t2-t1+1;
tt = (t1:t2)'.*delta_t - mid_t;

% % debug 
% tt = -2:2;
% d = [-2.5, -1.1, 0.4 1.5];
% mid_t = 0;
% delta_t = 1;
% Nreb = size(tt, 2);

fd = (d+mid_t)./delta_t;
dindex = floor(fd) - t1 + 2;
dindex(dindex<=0) = 1;
dindex(dindex>Nreb) = Nreb+1;
tindex = nan(Nreb+1, 1);
tindex(dindex) = 1:Npixel;
tindex = fillmissing(tindex(1:end-1), 'previous');

interalpha = (tt - d(tindex))./(d(tindex+1)-d(tindex));

% B = zeros(Nreb, Nview);
B = A(tindex,:).*repmat(1-interalpha, 1, Nview) + A(tindex+1,:).*repmat(interalpha, 1, Nview);

% B(700, 1) = 10000;




