function Aoff = offfocalconvHelical(A, detector, focalposition, Nviewprot, offwidth, offintensity, offedge)
% off-focal convolution (only for helical)
% Aoff = offfocalconvHelical(A, detector, focalposition, Nviewprot, offwidth, offintensity);
% A is the raw data on fan beam (one rotation)
% Aoff is the offfocal fix

SID = detector.SID;
SDD = detector.SDD;
DonL = SID/SDD;

% I know
[Npixel, Nslice, Nview] = size(A);
Nps = Npixel*Nslice;

% fan angles
y = detector.position(1:Npixel, 2) - focalposition(2);
x = detector.position(1:Npixel, 1) - focalposition(1);
fanangles = atan2(y, x);

% off-focal geodetic line
alpha = acos(DonL);
phi = fanangles-pi/2;
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);
% y1 = phi-log((cos(alpha)-cos(phi)-1i.*sin(phi).*sin(alpha))./(cos(alpha)*cos(phi)-1))./sqrt(cos(alpha)^2-1);
% I know phi_off = real(y1)

% d_view
delta_view = pi*2/Nviewprot;
% start-end view of off-focal
offstartview = 1 + floor(double(max(phi) - max(phi_off))/delta_view);
offendview = Nview + ceil(double(min(phi) - min(phi_off))/delta_view);
% Nviewoff = offendview - offstartview + 1;

% off-focal tau-measure
t_off = tan(phi - phi_off);
offsample = max(2^ceil(log2(Npixel)), 64);  % =1024
t = linspace(min(t_off), max(t_off), offsample);
% prepare the interpolation on the measure
[t_index, t_alpha] = interpprepare(t_off, t(:), 'extrap');

% off-focal kernel (sinc-window)
offwidth_nrm = offwidth/SID/(max(t_off)-min(t_off));
offkernel = offfocalsinckernel(offintensity, offwidth_nrm, offedge, offsample);

% ini Aoff
% Aoff1 = zeros(Npixel, Nviewoff);

% to interplate the raw data to off-focal geodetic line
f = (phi_off - phi)./delta_view;
intp_idx0 = double(floor(f));
intp_alpha = repmat(f-intp_idx0, Nslice, 1);
% off-focal view index
index_v = offstartview : offendview;
% neighboring data view index
intp_idx1 = intp_idx0 + index_v - 1;
intp_idx2 = intp_idx1 + 1;
intp_idx1(intp_idx1<0) = 0;
intp_idx1(intp_idx1>Nview-1) = Nview-1;
intp_idx2(intp_idx2<0) = 0;
intp_idx2(intp_idx2>Nview-1) = Nview-1;
% view index + pixiel index
intp_idx1 = repmat(intp_idx1.*Nps, Nslice, 1) + (1:Nps)';
intp_idx2 = repmat(intp_idx2.*Nps, Nslice, 1) + (1:Nps)';
% interp
Aoff1 = A(intp_idx1).*(1-intp_alpha) + A(intp_idx2).*intp_alpha;
Aoff1 = reshape(Aoff1, Npixel, []);

% to interplate the Aoff to tau-measure
Aoff1 = Aoff1(t_index(:,1),:).*t_alpha(:,1) + Aoff1(t_index(:,2),:).*t_alpha(:,2);

% convolution by off-focal kernel
Aoff1 = ifft(fft(Aoff1, offsample).*offkernel, 'symmetric');

% to interplate back the Aoff from tau-measure to off-focal geodetic line
[tb_index, tb_alpha] = interpprepare(t(:), t_off, 'extrap');
Aoff1 = Aoff1(tb_index(:,1),:).*tb_alpha(:,1) + Aoff1(tb_index(:,2),:).*tb_alpha(:,2);
Aoff1 = reshape(Aoff1, Nps, []);

% to interplate back to fan beam
Aoff = zeros(Npixel*Nslice, Nview);
Aoff(intp_idx1) = Aoff(intp_idx1) + Aoff1.*(1-intp_alpha);
Aoff(intp_idx2) = Aoff(intp_idx2) + Aoff1.*intp_alpha;

Aoff = reshape(Aoff, Npixel, Nslice, Nview);

end