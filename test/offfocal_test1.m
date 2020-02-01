focalposition = SYS.source.focalposition;
detector = SYS.detector;
SID = SYS.detector.detector_corr.SID;
SDD = SYS.detector.detector_corr.SDD;
DonL = SID/SDD;
Nviewprot = SYS.protocol.viewperrot;
A = Data.P{1};

tic;

offfocal_width = 50;
offfocal_intensity = 0.05;

Npixel = double(detector.Npixel);
Nslice = double(detector.Nslice);
% mid_U = single(detector.mid_U);
Nps = Npixel*Nslice;

% fan angles
y = detector.position(1:Npixel, 2) - focalposition(2);
x = detector.position(1:Npixel, 1) - focalposition(1);
fanangles = atan2(y, x);

alpha = acos(DonL);
phi = fanangles-pi/2;

% alpha = acos(0.5);

y1 = phi-log((cos(alpha)-cos(phi)-1i.*sin(phi).*sin(alpha))./(cos(alpha)*cos(phi)-1))./sqrt(cos(alpha)^2-1);
y2 = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);

t_off = tan(phi - y2);
offsample = max(2^ceil(log2(Npixel)), 64);
t = linspace(min(t_off), max(t_off), offsample);
[t_index, t_alpha] = interpprepare(t_off, t(:), 'extrap');
p_off = offfocal_width/SID/(max(t_off)-min(t_off));
tt = [0:offsample/2, -offsample/2+1:-1]';
offkernel = sinc(p_off.*tt).*offfocal_intensity;


delta_view = pi*2/Nviewprot;
f = (y2 - phi)./delta_view;
intp_idx0 = floor(f);
intp_a = f-intp_idx0;

index_v = 1:Nviewprot;
intp_idx1 = mod(intp_idx0 + index_v - 1, Nviewprot);
intp_idx2 = mod(intp_idx1, Nviewprot);

intp_idx1 = repmat(intp_idx1.*Nps, Nslice, 1) + (1:Nps)';
intp_idx2 = repmat(intp_idx2.*Nps, Nslice, 1) + (1:Nps)';
intp_a = repmat(intp_a, Nslice, 1);
A1 = A(intp_idx1).*(1-intp_a) + A(intp_idx2).*intp_a;
A1 = reshape(A1, Npixel, []);

% A2 = zeros(offsample, Nslice*Nviewprot);
A2 = A1(t_index(:,1),:).*t_alpha(:,1) + A1(t_index(:,2),:).*t_alpha(:,2);

A3 = ifft(fft(A2, offsample).*offkernel,'symmetric');
% A3 = A3(1:Npixel,:);

[tb_index, tb_alpha] = interpprepare(t(:), t_off, 'extrap');
A4 = A3(tb_index(:,1),:).*tb_alpha(:,1) + A3(tb_index(:,2),:).*tb_alpha(:,2);
A4 = reshape(A4, Nps, Nviewprot);

A5 = zeros(Nps, Nviewprot);
A5(intp_idx1) = A5(intp_idx1) + A4.*(1-intp_a);
A5(intp_idx2) = A5(intp_idx2) + A4.*intp_a;

toc;