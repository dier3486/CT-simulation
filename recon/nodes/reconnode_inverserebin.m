function [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status)
% cali node, inverse of the rebin
% [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status);
% use to find out the inverse of the rebin

% parameters to use
Nshot = prmflow.recon.Nshot;
Nviewprot = prmflow.recon.Nviewprot;
delta_view = prmflow.recon.delta_view;
focalspot = prmflow.system.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
detector = prmflow.system.detector;
Npixel = double(detector.Npixel);
Nslice = double(detector.Nmergedslice);
% mid_U = single(detector.mid_U);
Nreb = prmflow.rebin.Nreb;
delta_d = prmflow.rebin.delta_d;
midchannel = prmflow.rebin.midchannel;
% I know it is the midchannel after rebin
startvindex = prmflow.rebin.startvindex;
% if isfield(prmflow.rebin, 'Nleft')
%     Nleft = prmflow.rebin.Nleft;
% else
%     Nleft = 0;
% end

% NO QDO!

% fan angles
y = detector.position(1:Npixel, 2) - focalposition(2);
x = detector.position(1:Npixel, 1) - focalposition(1);
fanangles = atan2(y, x);

% invers rebin samples on theta-d space
f = fanangles./delta_view;
fv = mod(f+(0:Nviewprot-1), Nviewprot)+1;
dv = -detector.SID.*cos(fanangles)./delta_d + midchannel;
dv = repmat(dv, 1, Nviewprot);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, Nviewprot, Nshot);
% ini
D = zeros(Npixel, Nslice, Nviewprot, Nshot);
for ishot = 1:Nshot
    for islice = 1:Nslice
        A = squeeze(dataflow.rawdata(:, islice, :, ishot));
        A = [A(:, end-startvindex+2:end) A(:, 1:end-startvindex+1)];
        A = [A A(:,1)];
        D(:, islice, :, ishot) = interp2(A,fv,dv, 'linear', 0);
    end
end

% to return
dataflow.rawdata = D;
prmflow.recon.Npixel = Npixel;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end