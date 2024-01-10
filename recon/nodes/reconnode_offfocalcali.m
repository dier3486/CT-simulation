function [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status)
% recon node, off-focal calibration based on scanning a block column 
% [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status);
% For perfect off-focal correction.

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% parameters set in pipe
caliprm = prmflow.pipe.(status.nodename);

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
    % the corrversion shall >= v2.0, the version 1.x are not supported to do the calibration base on data.
else
    corrversion = 'v2.0';
end

Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
Nviewprot = prmflow.raw.Nviewprot;
Nview = size(dataflow.rawdata, 2);

if isfield(caliprm, 'offsample')
    Noffsample = caliprm.offsample;
    % suggest to set the Noffsample in 2048
else
    Noffsample = max(2^ceil(log2(Npixel)), 64);  % =1024
end
if isfield(caliprm, 'tauremeasure')
    at3_scale = caliprm.tauremeasure.t3;
else
    at3_scale = 0.1;
end

% special data correction
% I know the air calibration table is
aircorr = prmflow.corrtable.Air;
airmean = mean(reshape(aircorr.main, [], aircorr.Nsection), 2);
minttime = mean(dataflow.rawhead.Integration_Time);
% inverse the log2 data to intensity
airfix = 2.^(-airmean + log2(minttime));
raw0 = 2.^(-dataflow.rawdata - airmean + log2(minttime));
airfix = reshape(airfix, Npixel, Nslice);
raw0 = reshape(raw0, Npixel, Nslice, Nview);

% dirty correction for low-signal issue
% L-H correction
raw0 = (raw0(:, 1:2:end, :) - raw0(:, 2:2:end, :))./(airfix(:, 1:2:end)-airfix(:, 2:2:end));

% focalposition
focalspot = prmflow.raw.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% no DFS!

% detector
detector = prmflow.system.detector;

% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice);

% half slice
Nslice2 = Nslice/2;
fanangles2 = (fanangles(:, 1:2:end) + fanangles(:, 2:2:end))./2;


SID = detector.SID;
SDD = detector.SDD;
alpha = acos(SID/SDD);
phi2 = fanangles2-pi/2;
phi2_off = phi2 - atan( sin(phi2).*sin(alpha)./(cos(phi2)-cos(alpha)) ) ./ sin(alpha);

t_off = double(phi2 - phi2_off);
dt = (max(t_off(:)) - min(t_off(:)))/(Noffsample - 1);
t0 = linspace(min(t_off(:)), max(t_off(:)), Noffsample)';

% tau-measure scale
tb = max(abs(t_off(:)));
t_resp = (t0 + t0.^3.*at3_scale).*(tb/(tb + tb^3*at3_scale));
dt_resp = (1 + t0.^2.*at3_scale.*3).*(tb/(tb + tb^3*at3_scale));
% dt_resp_off = (1 + t_off.^2.*at3_scale.*3).*(tb/(tb + tb^3*at3_scale));
t_odd = t0;

index_p = single(1:Npixel)';
index_v = 1:Nviewprot;
delta_view = pi*2/Nviewprot;
Dphi2 = phi2 - phi2_off;

% interp raw0 to off-focal measure space
raw1 = zeros(Noffsample, Nslice2, Nviewprot);
for islice = 1:Nslice2
    Df = mod(index_v - Dphi2(:, islice)./delta_view - 1, Nviewprot) + 1;
    raw0_isl = [squeeze(raw0(:, islice, :)) raw0(:, islice, 1)];
    rawtmp = interp2(raw0_isl , Df, repmat(index_p, 1, Nviewprot));

    raw1(:, islice, :) = interp1(t_off(:, islice), rawtmp, t_resp(:), 'linear', 'extrap');
end
% mean (selected slices) of raw1 to do the curve fitting
raw2 = squeeze(mean(raw1(:, 2:end-1, :), 2));
% skip the edge slices

% prepare for curve fitting data sample based on raw2
data0 = ones(Noffsample, Nviewprot);
data_d = zeros(Noffsample, Nviewprot);
Wd = ones(Noffsample, Nviewprot);
% to find out the correction target (data0) and the error to fix (data_d)
if isfield(caliprm, 'edgesample')
    edgesample = caliprm.edgesample;
else
    edgesample = 40;
end
mk = round(edgesample/2);

options = optimoptions('lsqnonlin','Display','off');
% adpative weighting
if isfield(caliprm, 'weightingscale')
    Wscale = caliprm.weightingscale;
else
    Wscale = 12.0;
end
U = zeros(4, Nviewprot);
for iview = 1:Nviewprot
    index_l = find(raw2(:, iview)<0.5, 1, 'first');
    index_r = find(raw2(:, iview)<0.5, 1, 'last');
    indexLR = (max(index_l-edgesample,1) : min(index_r+edgesample, Noffsample))';
    data_r = raw2(indexLR, iview);
    umin = mean(mink(data_r, mk));
    u = lsqnonlin(@(u) offwell(indexLR, u, umin) - data_r, [index_l, index_r, 1, 1], [], [], options);
    U(:, iview) = u(:);
    ufix = u;
%     ufix([3 4]) = edgeshapen;
%     Ru = offwell(indexLR, u, umin);  % use to observe the curves
    Rfix = offwell(indexLR, ufix, umin);
    data0(indexLR, iview) = Rfix;
    data_d(indexLR, iview) = data_r - Rfix;
    
    Gw = exp(-(indexLR-ufix(1)).^2.*ufix(3)^2).*(ufix(3)/sqrt(pi)) + exp(-(indexLR-ufix(2)).^2.*ufix(4)^2).*(ufix(4)/sqrt(pi));
    Wd(indexLR, iview) = exp(-Gw.*Wscale);
end
% the data0 and data_d shall be double, no single float plz.
1; % block here to observe the correction target and raw data to verify those parameters

if isfield(caliprm, 'edgeshapen')
    edgeshapen = caliprm.edgeshapen;
else
    edgeshapen = mean(mean(U(3:4, :)));
end

% rescale the data_d by dtau
data_d = data_d.*dt_resp;


Fd = fft(data_d);
Fd(Noffsample/2+1, :) = real(Fd(Noffsample/2, :));
% Fa = fft(data2);
% Fa(Noffsample/2+1, :) = real(Fa(Noffsample/2, :));
Fa0 = fft(data0);
Fa0(Noffsample/2+1, :) = real(Fa0(Noffsample/2, :));

Fw = fft(Wd);

% A, B shall in double
A = zeros(Noffsample, Noffsample);
index_k = [(0 : Noffsample/2)  (-Noffsample/2+1) : -1];
B = zeros(Noffsample, 1);
index_w = mod(index_k - index_k(:), Noffsample) + 1;

for ii = 1:Nviewprot
    W_ii = reshape(Fw(index_w(:), ii), Noffsample, Noffsample);
    WF_ii = W_ii.*Fa0(:, ii).';
    A = A + WF_ii'*WF_ii;
    B = B + WF_ii'*(W_ii*Fd(:, ii));
end

% raw kernel A\B

% % normal style
% L = diag(index_k.^2.*(pi*2/Noffsample)^2.*(sqrt(Nviewprot)/2));
% L(Noffsample/2+1, Noffsample/2+1) = L(Noffsample/2+1, Noffsample/2+1)./2;
% fk0 = A\B;
% k0 = real(ifft(fk0));

% symatric style (need to simplify codes later)
A1 = [A(:, 1) (A(:,2:Noffsample/2) + conj(A(:,Noffsample:-1:Noffsample/2+2))) A(:, Noffsample/2+1)];    
A1 = [A1(1, :); (A1(2:Noffsample/2, :) + conj(A1(Noffsample:-1:Noffsample/2+2, :)))./2; A1(Noffsample/2+1, :)];
B1 = [B(1); (B(2:Noffsample/2) + conj(B(Noffsample:-1:Noffsample/2+2)))./2; B(Noffsample/2+1)];

L2 = diag(index_k(1:Noffsample/2+1).^2.*(pi*2/Noffsample)^2.*(sqrt(Nviewprot)/2)).*2;
L2(end, end) = L2(end, end)./4;

if isfield(caliprm, 'Laplacescale')
    Lscale = caliprm.Laplacescale;
else
    Lscale = 1.0;
end
% the effect of L is very small, but can clear the singularity

% solve (A+L)*fk = B
fk2 = (real(A1) + L2.*Lscale)\real(B1);

G0 = exp(-index_k(:).^2.*edgeshapen.^2).*(edgeshapen/sqrt(pi));
FG0 = exp(-(index_k(:)./Noffsample*2*pi).^2.*(1/edgeshapen^2/4));
rawkernel = ifft(([fk2; fk2(Noffsample/2:-1:2)]+1).*FG0, 'symmetric');
1; % block here to observe the off-focal kernel

mk = Noffsample/32;
g2 = rawkernel(1:mk+1);

% try to smooth the center reagion
[~, smax] = max(g2(4:mk+1) - G0(4:mk+1));
smax = max(smax + 3 - 2, 2);
x_s = (smax-1:mk)';
u2 = lsqnonlin(@(u) ((exp(-x_s.^2./2./u(1)^2).*u(2)).^2 - g2(smax:end).^2).*Noffsample^2, [smax g2(smax)]);
Gu2 = exp(-(0:mk)'.^2./2./u2(1)^2).*u2(2);
Wg = exp(-G0(1:mk+1).*Noffsample.*2.0);
g2 = Gu2.*(1-Wg) + g2.*Wg;
g2(g2<0) = 0;

if isfield(caliprm, 'Iscale')
    Iscale = caliprm.Iscale;
else
    Iscale = 1.0;
end
if isfield(caliprm, 'Iodd')
    Iodd = caliprm.Iodd;
else
    Iodd = 0.6;
end

% % to observe the correction
% g2s = zeros(Noffsample, 1);
% g2s(1:mk+1) = g2 + g2.*1i;
% g2s(end-mk+1:end) = flipud(g2(2:end) - g2(2:end).*1i);

% Koff0 = fft(g2s);
% Koff1 = Koff0 - Koff0(1);
% % the Taylor
% Koff2 = Koff1+Koff1.^2;
% Koff3 = Koff1+Koff1.^2+Koff1.^3;
% Koff_inf = 1./(1-Koff1)-1; 
% 
% % select one of them
% Koff_work = Koff_inf;
% 
% % fix on raw0
% Dfix0 = ifft(fft(reshape(raw1, Noffsample, [])).*Koff_work)./dt_resp;
% Dfix0 = reshape(Dfix0, Noffsample, Nslice2, []);
% 
% Dfix0_m0 = zeros(Npixel, Nslice2, Nviewprot);
% for islice = 1:Nslice2
%     Dfb = mod(index_v + Dphi2(:, islice)./delta_view - 1, Nviewprot) + 1;
%     Dfix0_tmp = interp1(t_resp, squeeze(Dfix0(:, islice, :)), Dphi2(:, islice), 'linear', 0);
%     Dfix0_m0(:, islice, :) = interp2([Dfix0_tmp Dfix0_tmp(:,1)] , Dfb, repmat(index_p, 1, Nviewprot));
% end
% Dfix0_m0 = Dfix0_m0.*Iscale;
% 
% taodd2 = interp1(t0, t_odd, Dphi2, 'linear','extrap').*Iodd;
% raw0_fix = raw0 - real(Dfix0_m0) - imag(Dfix0_m0).*taodd2;
1; % block here to observe the correction
% plz try the calibation parameters to make the raw0_fix best.

if isfield(caliprm, 'crossrate')
    crossrate = caliprm.crossrate;
else
    crossrate = 0;
end

% offfocalcorr
offfocalcorr = caliprmforcorr(prmflow, corrversion);

offkernel = struct();
offkernel.edgesample = edgesample;
offkernel.Noffsample = Noffsample;
offkernel.weightingscale = Wscale;
offkernel.Laplacescale = Lscale;
offkernel.tauremeasure_at3 = at3_scale;
offkernel.t0 = t0;
offkernel.tresp = t_resp;
offkernel.tscale = 1./dt_resp;
offkernel.tscaleodd = t_odd;
offkernel.edgeshapen = edgeshapen;

offkernel.x = (-mk:mk)'.*(SID*dt);
offkernel.dx = SID*dt;
offkernel.curve = [flipud(g2(2:end)); g2];
offkernel.curveodd = [-flipud(g2(2:end)); 0; g2(2:end)];
offkernel.Iscale = Iscale;
offkernel.Iodd = Iodd;
offkernel.crossrate = crossrate;

% merge
offfocalcorr = structmerge(offfocalcorr, offkernel);

% steps
if isfield(offfocalcorr, 'offintensity')
    offfocalcorr.stepnumber = length(offfocalcorr.offintensity);
else
    offfocalcorr.stepnumber = 1;
end

% to return 
dataflow.offfocalcorr =  offfocalcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

function y = offwell(x, u, umin)
y = ((1-erf((x-u(1)).*u(3))./2 - erf((-x+u(2)).*u(4))./2)).*(1-umin) + umin;

end