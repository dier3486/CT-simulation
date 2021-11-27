% data
load('D:\matlab\Data\SINO\offfocal\offin.mat');

Nview = prmflow.recon.Nview;
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nviewprot = prmflow.recon.Nviewprot;
scantype = prmflow.recon.scan;

% fix to QFS
focalposition = prmflow.system.focalposition(1, :);
% detector
detector = prmflow.system.detector;
SID = detector.SID;
SDD = detector.SDD;
% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
% mean
fanangles = mean(reshape(fanangles, Npixel, Nslice), 2);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% air rate
% dataflow.rawdata = dataflow.rawdata+offcorr.airrate;
% exp
dataflow.rawdata = 2.^(-dataflow.rawdata);
% slice mean
dataflow.rawdata = squeeze(mean(reshape(dataflow.rawdata, Npixel, Nslice, Nview), 2));

crossrate = 0;
% prepare the off-focal base intensity with a z-crossed method
% Aoff = offfocalzcross(reshape(dataflow.rawdata, Npixel, Nslice, Nview), crossrate);
Aoff = reshape(dataflow.rawdata, Npixel, Nview);

% fill 0
Aoff = fillmissing(Aoff, 'constant',0);

% A0
A0 = (Aoff>0.5).*1.0;

%--% off conv

% reshape
Aoff = reshape(Aoff, Npixel, Nviewprot);
Np = Npixel;

% off-focal geodetic line
DonL = SID/SDD;
alpha = acos(DonL);
phi = fanangles-pi/2;
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);

% off-focal tau-measure
t_off = tan(phi - phi_off);
offsample = max(2^ceil(log2(Npixel)), 64);  % =1024
t = linspace(min(t_off), max(t_off), offsample);
% prepare the interpolation on the measure
[t_index, t_alpha] = interpprepare(t_off, t(:), 'extrap');

% to interplate the raw data to off-focal geodetic line
delta_view = pi*2/Nviewprot;
f = (phi_off - phi)./delta_view;
intp_idx0 = double(floor(f));
intp_alpha = f-intp_idx0;
% intp_alpha = repmat(intp_alpha, Nslice, 1);
% to interplate on all the views
index_v = 1:Nviewprot;
intp_idx1 = mod(intp_idx0 + index_v - 1, Nviewprot);
intp_idx2 = mod(intp_idx1 + 1, Nviewprot);
intp_idx1 = intp_idx1.*Np + (1:Np)';
intp_idx2 = intp_idx2.*Np + (1:Np)';
Aoff1 = Aoff(intp_idx1).*(1-intp_alpha) + Aoff(intp_idx2).*intp_alpha;
Aoff1 = reshape(Aoff1, Npixel, []);
Aoff0 = A0(intp_idx1).*(1-intp_alpha) + A0(intp_idx2).*intp_alpha;
Aoff0 = reshape(Aoff0, Npixel, []);

% to interplate the Aoff to tau-measure
Aoff1 = Aoff1(t_index(:,1),:).*t_alpha(:,1) + Aoff1(t_index(:,2),:).*t_alpha(:,2);
Aoff0 = Aoff0(t_index(:,1),:).*t_alpha(:,1) + Aoff0(t_index(:,2),:).*t_alpha(:,2);

% W
Nv = size(Aoff0, 2);
W = 1 - conv2([zeros(1, Nv); diff(Aoff0>0)~=0], ones(6,1), 'same');

% K
K = [0:offsample/2  -offsample/2+1:-1]'.*(2*pi/offsample);
% fft
Boff1 = fft(Aoff1, offsample);
Boff0 = fft(Aoff0, offsample);

% 1
b1 = mean((Boff0-Boff1).*conj(Boff1), 2)./offsample;
alpha = 2.0;
A1 = mean(conj(Boff1).*Boff1, 2)./offsample^2 + K.^2.*alpha./offsample^2;

y1 = ifft(b1./A1);
y1 = [y1(offsample/2+2:end); y1(1:offsample/2+1)];

% 2
Bv = fft((Aoff0-Aoff1).*W, offsample);
Wfft = fft(W, offsample);

N1 = offsample/2+1;
b2 = mean(Bv.*conj(Boff1), 2)./offsample;

% w_i-j
windex1 = (-offsample/2+1:offsample/2)' - (-offsample/2+1:offsample/2);
windex1(abs(windex1)>=offsample/2) = offsample/2;
windex1 = mod(windex1, offsample) + 1;

% % w_i-j
% windex1 = (0:offsample/2)' + (0:-1:-offsample/2);
% windex1 = mod(windex1, offsample) + 1;
% % w_i+j
% windex2 = (0:offsample/2)' + (0:1:offsample/2);
% windex2(windex2>offsample/2) = offsample/2;
% windex2 = windex2 + 1;
% % I know Wfft(offsample/2+1)=0.
% 
% M = zeros(N1, N1);
% for ii = 1:Nv
%     W_ii = Wfft(:, ii);
%     M1 = (Boff1(1:N1, ii)*Boff1(1:N1, ii)').*conj(W_ii(windex1));
%     M2 = (Boff1(1:N1, ii)*Boff1(1:N1, ii).').*conj(W_ii(windex2));
%     M = M + (real(M1)+real(M2));
% end
% M = M./Nv;
