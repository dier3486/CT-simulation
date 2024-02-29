function Aoff = offfocalconvAxial(A, fanangles, SID, SDD, Nviewprot, offwidth, offintensity, offedge)
% off-focal convolution (only for axial)
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
phi = fanangles-pi/2;
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);
% y1 = phi-log((cos(alpha)-cos(phi)-1i.*sin(phi).*sin(alpha))./(cos(alpha)*cos(phi)-1))./sqrt(cos(alpha)^2-1);
% I know phi_off = real(y1)

% off-focal tau-measure
t_off = tan(phi - phi_off);
offsample = max(2^ceil(log2(Npixel)), 64);  % =1024
t = linspace(min(t_off), max(t_off), offsample);
% prepare the interpolation on the measure
[t_index, t_alpha] = interpprepare(t_off, t(:), 'extrap');

% off-focal kernel (sinc-window)
offwidth_nrm = offwidth/SID/(max(t_off)-min(t_off));
offkernel = offfocalsinckernel(offintensity, offwidth_nrm, offedge, offsample);
% p_off = offwidth/SID/(max(t_off)-min(t_off));
% tt = [0:offsample/2, -offsample/2+1:-1]';
% offkernel = sinc(p_off.*tt).*offintensity;
% 
% % edge smooth
% edge_smooth = 1./(1+exp((-offedge+abs(tt)./(offsample/2)).*(10/(0.5-abs(offedge-0.5)))));
% edge_smooth = fillmissing(edge_smooth, 'nearest');
% offkernel = offkernel.*edge_smooth;

% to interplate the raw data to off-focal geodetic line
if Nviewprot>1
    delta_view = pi*2/Nviewprot;
    f = (phi_off - phi)./delta_view;
    intp_idx0 = double(floor(f));
    intp_alpha = f-intp_idx0;
    intp_alpha = repmat(intp_alpha, NAslice, 1);
    % to interplate on all the views
    index_v = 1:Nviewprot;
    intp_idx1 = mod(intp_idx0 + index_v - 1, Nviewprot);
    intp_idx2 = mod(intp_idx1 + 1, Nviewprot);
    intp_idx1 = repmat(intp_idx1.*Npa, NAslice, 1) + (1:Npa)';
    intp_idx2 = repmat(intp_idx2.*Npa, NAslice, 1) + (1:Npa)';
    Aoff1 = A(intp_idx1).*(1-intp_alpha) + A(intp_idx2).*intp_alpha;
    Aoff1 = reshape(Aoff1, Npixel, []);
else
    % Nviewprot==1 (for air simulation, NOT topo!)
    Aoff1 = reshape(A, Npixel, []);
end

% to interplate the Aoff to tau-measure
Aoff1 = Aoff1(t_index(:,1),:).*t_alpha(:,1) + Aoff1(t_index(:,2),:).*t_alpha(:,2);

% convolution by off-focal kernel
Aoff1 = ifft(fft(Aoff1, offsample).*offkernel,'symmetric');

% to interplate back the Aoff from tau-measure to off-focal geodetic line
[tb_index, tb_alpha] = interpprepare(t(:), t_off, 'extrap');
Aoff1 = Aoff1(tb_index(:,1),:).*tb_alpha(:,1) + Aoff1(tb_index(:,2),:).*tb_alpha(:,2);
Aoff1 = reshape(Aoff1, Npa, Nviewprot);

% to interplate back to fan beam
if Nviewprot>1
    Aoff = zeros(Npa, Nviewprot);
    Aoff(intp_idx1) = Aoff(intp_idx1) + Aoff1.*(1-intp_alpha);
    Aoff(intp_idx2) = Aoff(intp_idx2) + Aoff1.*intp_alpha;
else
    % Nviewprot==1 (for air simulation)
    Aoff = Aoff1;
end

% reshape as the input
Aoff = reshape(permute(reshape(Aoff, Npixel*Nslice, Nshot, Nviewprot), [1 3 2]), Npixel, Nslice, Nview);

end