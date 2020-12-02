% 'slope' rebin test code
% Atilt transform

% load('E:\data\simulation\TM\test\tiltrb_test1.mat');

Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);

gantrytilt = prmflow.protocol.gantrytilt*(pi/180);

Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;

% detector
detector = prmflow.system.detector;

delta_d = detector.hx_ISO;
dfan = atan(detector.hx_ISO/detector.SID);
% dfan = detector.hx_ISO/detector.SID;
idealfan = dfan.*((1:Npixel)-detector.mid_U);
t0 = 1-mod(detector.mid_U, 1);
t1 = ceil(tan(min(idealfan))*detector.SID/delta_d - t0);
t2 = floor(tan(max(idealfan))*detector.SID/delta_d - t0);
idealphi = atan(((t1:t2)+t0)./detector.SID);

Zdet = mean(reshape(detector.position(:,3), Npixel, Nslice));
Yshift = repelem(-Zdet(:).*tan(gantrytilt), Npixel, 1);

viewangle = dataflow.rawhead.viewangle;

dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

for iview = 1:Nview
    Rv = [cos(viewangle(iview))  sin(viewangle(iview));
         -sin(viewangle(iview))  cos(viewangle(iview))];
    focal_iview = focalposition(:, 1:2)*Rv;
    detpos_iview = detector.position(:, 1:2)*Rv;
    detpos_iview(:,2) = detpos_iview(:, 2) + Yshift;
    
    [fanangles, focalangle] = detpos2fanangles(detpos_iview, focal_iview);
    fanangles = reshape(fanangles, Npixel, Nslice);
    for islice = 1:Nslice
    end
end

% idealfan = mean(fanangles-pi/2)./(Npixel/2+1/2-detector.mid_U).*((1:Npixel)-detector.mid_U);


% % rebin prepare
% [prmflow, ~] = reconnode_rebinprepare(prmflow, status);
% 
% % Azi rebin
% [dataflow, prmflow, ~] = reconnode_Azirebin(dataflow, prmflow, status);