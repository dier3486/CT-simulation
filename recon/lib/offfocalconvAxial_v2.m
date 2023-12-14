function Aoff = offfocalconvAxial_v2(A, fanangles, SID, SDD, Nviewprot)
% off-focal convolution (only for axial) version2, in working
% Aoff = offfocalconvAxial(A, fanangles, SID, SDD, Nviewprot, offwidth, offintensity);
% A is the raw data on fan beam (one rotation)
% Aoff is the offfocal fix

% SID = detector.SID;
% SDD = detector.SDD;
DonL = SID/SDD;

% I know
[Npixel, Nslice, Nview] = size(A);
Nshot = Nview/Nviewprot;

% to support multi shots
A = reshape(A, Npixel*Nslice, Nviewprot, []);
A = reshape(permute(A, [1 3 2]), [], Nviewprot);
Npa = size(A, 1);
% I know Npa = Npixel*Nslice*Nshot;
NAslice = Nslice*Nshot;

% % fan angles
% y = detector.position(1:Npixel, 2) - focalposition(2);
% x = detector.position(1:Npixel, 1) - focalposition(1);
% fanangles = atan2(y, x);

% off-focal geodetic line
alpha = acos(DonL);
phi = fanangles(:)-pi/2;
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);
% y1 = phi-log((cos(alpha)-cos(phi)-1i.*sin(phi).*sin(alpha))./(cos(alpha)*cos(phi)-1))./sqrt(cos(alpha)^2-1);
% I know phi_off = real(y1)

% load off-focal kernel and measure
offkernel = load('F:\咨询\WY\3.0\offfocal所需数据\offfocal所需数据\低能\offkernel.mat');

% off-focal tau-measure
offsample = max(2^ceil(log2(Npixel)), 64);  % =1024
t_off = phi - phi_off;

% re-scaled measure
t0 = linspace(min(t_off), max(t_off), offsample)';
t_resp = interp1(offkernel.t0(:), offkernel.t_resp(:), t0, 'linear', 'extrap');
Krescale = interp1(offkernel.t0(:), offkernel.rescale(:), t_off, 'linear', 'extrap');
Krescale_odd = interp1(offkernel.t0(:), offkernel.rescale_odd(:), t_off, 'linear', 'extrap');

% off-focal kernel
dt0 = (max(t_off) - min(t_off))/(offsample - 1)*SID;
m = ceil(max(offkernel.x)/dt0);

kernelCurve = spectrumresample(offkernel.x, offkernel.curve, (-m:m).*dt0);
kernelK = zeros(offsample, 1);
kernelK(1:m+1) = kernelCurve(m+1:end);  kernelK(end-m+1:end) = kernelCurve(1:m);
kernelK = 1 - 1./(fft(kernelK) + 1);
kernelK = kernelK - kernelK(1);

% prepare for interp
index_p = single(1:Npixel)';
index_v = 1:Nviewprot;
delta_view = pi*2/Nviewprot;
Dphi = t_off./delta_view;

t1 = interp1(t_off, index_p, t_resp, 'linear', 'extrap');

% to interplate the raw data to off-focal geodetic line
if Nviewprot>1
    Df = mod(index_v - Dphi - 1, Nviewprot) + 1;
    Aoff1 = interp2([A A(:,1)] , Df, repmat(index_p, 1, Nviewprot));
else
    % Nviewprot==1 (for air simulation, NOT topo!)
    Aoff1 = reshape(A, Npixel, []);
end
% to interplate the Aoff to tau-measure
Aoff1 = interp1(Aoff1, t1, 'linear', 'extrap');

% convolution by off-focal kernel
Aoff1 = ifft(fft(Aoff1).*kernelK);

% to interplate back the Aoff from tau-measure to off-focal geodetic line
Aoff1 = interp1(t_resp, Aoff1, t_off, 'linear', 0);

% to interplate back to fan beam
if Nviewprot>1
    Dfb = mod(index_v + Dphi - 1, Nviewprot) + 1;
    Aoff = interp2([Aoff1 Aoff1(:,1)] , Dfb, repmat(index_p, 1, Nviewprot));
else
    % Nviewprot==1 (for air simulation)
    Aoff = Aoff1;
end

% rescale
Aoff = Aoff.*Krescale;
Aoff = real(Aoff) + imag(Aoff).*Krescale_odd;

% reshape as the input
Aoff = reshape(permute(reshape(Aoff, Npixel*Nslice, Nshot, Nviewprot), [1 3 2]), Npixel, Nslice, Nview);

end