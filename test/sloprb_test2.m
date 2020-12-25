% 'slope' rebin test code
% Atilt transform

load('E:\data\simulation\TM\test\tiltrb_test1.mat');

% Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);

gantrytilt = prmflow.protocol.gantrytilt*(pi/180);

Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
% detector
detector = prmflow.system.detector;

delta_d = detector.hx_ISO;
dfan = atan(detector.hx_ISO/detector.SID);
% dfan = detector.hx_ISO/detector.SID;
idealfan = dfan.*((1:Npixel)'-detector.mid_U);
t0 = mod(detector.mid_U, 1);
t1 = ceil(sin(min(idealfan))*detector.SID/delta_d + t0);
t2 = floor(sin(max(idealfan))*detector.SID/delta_d + t0);
idealphi = asin(((t1:t2)'-t0).*(delta_d/detector.SID));
Nphi = t2-t1+1;
midU_phi = -t1+1+t0;

Zdet = mean(reshape(detector.position(:,3), Npixel, Nslice));
Yshift = repelem(-Zdet(:).*tan(gantrytilt), Npixel, 1);

viewangle = dataflow.rawhead.viewangle;

dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

tic;
data1 = zeros(Npixel, Nslice, Nview);
data2 = zeros(Nphi, Nslice, Nview);
for iview = 1:Nview
    Rv = [cos(viewangle(iview))  sin(viewangle(iview));
         -sin(viewangle(iview))  cos(viewangle(iview))];
    focal_iview = focalposition(:, 1:2)*Rv;
    detpos_iview = detector.position(:, 1:2)*Rv;
    detpos_iview(:,2) = detpos_iview(:, 2) + Yshift;
    
    [fanangles, focalangle] = detpos2fanangles(detpos_iview, focal_iview);
    fanangles = reshape(fanangles, Npixel, Nslice) - focalangle;
    fanangles = mod(fanangles+pi, pi*2) - pi;
    for islice = 1:Nslice
        data1(:, islice, iview) = interp1(fanangles(:, islice), dataflow.rawdata(:, islice, iview), idealfan, 'linear', 0);
        data2(:, islice, iview) = interp1(fanangles(:, islice), dataflow.rawdata(:, islice, iview), idealphi, 'linear', 0);
    end
%     if iview>=288
%         1;
%     end
    % z interp
    
end
toc;

% 1
prmflow.recon.fanangles = repmat(idealfan, 1, Nslice);
prmflow.recon.focalangle = 0;
dataflow.rawdata = data1;
[prmflow, ~] = reconnode_rebinprepare(prmflow, status);
[dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);
[dataflow, prmflow, ~] = reconnode_Radialrebin(dataflow, prmflow, status);
data1 = reshape(dataflow.rawdata, Nphi, Nslice, Nview);
rebin1 = prmflow.rebin;
recon1 = prmflow.recon;

% 2
prmflow.recon.fanangles = repmat(idealphi, 1, Nslice);
prmflow.recon.focalangle = 0;
prmflow.recon.Npixel = Nphi;
prmflow.system.detector.mid_U = midU_phi;
dataflow.rawdata = data2;
[prmflow, ~] = reconnode_rebinprepare(prmflow, status);
[dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);
data2 = reshape(dataflow.rawdata, Nphi, Nslice, Nview);
rebin2 = prmflow.rebin;
recon2 = prmflow.recon;
