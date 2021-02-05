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
Nshot = prmflow.recon.Nshot;
% detector
detector = prmflow.system.detector;

[fanangles0, focalangle0] = detpos2fanangles(detector.position, focalposition);
fanangles0 = reshape(fanangles0, Npixel, Nslice) - pi/2;

delta_d = detector.hx_ISO;
dfan = atan(double(detector.hx_ISO)/double(detector.SID));
% dfan = detector.hx_ISO/detector.SID;

[dfan1, dfan2] = idealdfan(double(mean(fanangles0, 2)), double(detector.mid_U));

% idealfan = dfan.*((1:Npixel)'-detector.mid_U);
idealfan = dfan1.*((1:Npixel)'-(detector.mid_U));
% idealfan = single(idealfan);
% idealfan = linspace(dfan2.*(1-detector.mid_U), dfan2.*(Npixel-detector.mid_U), Npixel)';
t0 = mod(detector.mid_U, 1);
t1 = ceil(sin(min(idealfan))*detector.SID/delta_d + t0);
t2 = floor(sin(max(idealfan))*detector.SID/delta_d + t0);
idealphi = asin(((t1:t2)'-t0).*(double(delta_d)/double(detector.SID)));
Nphi = t2-t1+1;
midU_phi = -t1+1+t0;

Zdet = mean(reshape(detector.position(:,3), Npixel, Nslice));
Yshift = repelem(-Zdet(:).*tan(gantrytilt), Npixel, 1);
Yshift1 = -Zdet.*tan(gantrytilt)./detector.SDD;

delta_z = detector.hz_ISO;
Zgrid0 = single(1:Nslice);
Zsec = cos(gantrytilt)./cos(idealphi);
Zsec(Zsec>1) = 1;
Zgrid = Zsec*(-(Nslice-1)/2 : (Nslice-1)/2) + (Nslice+1)/2;
Xgrid0 = 1:Nphi;
Xgrid = repmat(Xgrid0(:), 1, Nslice);



faninterpkern = zeros(Npixel, Nslice, 'single');
for islice = 1:Nslice
    faninterpkern(:, islice) = interp1(fanangles0(:,islice), 1:Npixel, idealfan, 'linear', 'extrap');
end

viewangle = dataflow.rawhead.viewangle;

dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

tic;
data0 = zeros(Npixel, Nslice, Nview);
data1 = zeros(Nphi, Nslice, Nview);
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
    
    phi1 = inverseyshift(idealphi, viewangle(iview), Yshift1);
    phi1 = reshape(phi1, Nphi, Nslice);
    phi1_index = phi1./dfan + detector.mid_U;
    for islice = 1:Nslice
% %         data1(:, islice, iview) = interp1(fanangles(:, islice), dataflow.rawdata(:, islice, iview), idealfan, 'linear', 0);
%         data1(:, islice, iview) = interp1(fanangles(:, islice), dataflow.rawdata(:, islice, iview), idealphi, 'linear', 0);
        data_tmp = interp1(dataflow.rawdata(:, islice, iview), faninterpkern(:, islice), 'linear', 'extrap');
        data2(:, islice, iview) = interp1(data_tmp, phi1_index(:, islice), 'linear', 0);
        data1(:, islice, iview) = interp1(fanangles0(:, islice), dataflow.rawdata(:, islice, iview), phi1(:, islice), 'linear', 0);
    end

%     if iview>=288
%         1;
%     end
    % z interp
    % ignore the components of Yshift;
    data1(:, :, iview) = interp2(Zgrid0, Xgrid0, data1(:, :, iview), Zgrid, Xgrid);
%     data2(:, :, iview) = interp2(Zgrid0, Xgrid0, data2(:, :, iview), Zgrid, Xgrid);
    data2(:, :, iview) = interp2(data2(:, :, iview), Zgrid, Xgrid);
end
toc;

% 2
prmflow.recon.fanangles = repmat(idealphi, 1, Nslice);
prmflow.recon.focalangle = 0;
prmflow.recon.Npixel = Nphi;
prmflow.system.detector.mid_U = midU_phi;
dataflow.rawdata = data2;
[prmflow, ~] = reconnode_rebinprepare(prmflow, status);
[dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);
% prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = prmflow.rebin.Nviewprot;
prmflow.recon.Nview = prmflow.rebin.Nviewprot*Nshot;
prmflow.recon.midchannel = prmflow.rebin.midchannel;
prmflow.recon.delta_d = prmflow.rebin.delta_d;
prmflow.recon.delta_z = prmflow.rebin.delta_z;
prmflow.recon.SID = prmflow.rebin.SID;

% Filter
status.nodename = 'Filter';
[dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);

% BP prepare
status.nodename = 'Backprojection';
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

function phi1 = inverseyshift(phi2, theta, Ys)

% a = tan(phi2).^2+1;
% b = -2.*(cos(theta).*tan(phi2) + sin(theta)).*Ys;
% c = (tan(phi2).*cos(theta) + sin(theta)).^2.*Ys.^2-tan(phi2).^2;
% 
% phi1 = asin((-b+sqrt(b.^2-a.*c.*4))./a./2);

phi1 = phi2 + asin(Ys.*sin(theta+phi2));

end

function [dfan1, dfan2] = idealdfan(fanangles, mid_U)

N = length(fanangles);
index = (1:N) - mid_U;

dfan1 = (index*fanangles(:)) / sum(index.^2);

f2 = fanangles(:) - index(:).*dfan1;
dfan2 = (index*f2) / sum(index.^2);

end
