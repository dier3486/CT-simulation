function [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status)
% recon node, interp the rawdata to equal-fanangle
% [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status);
%

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nfocal = prmflow.recon.Nfocal;
focalspot = prmflow.recon.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);

% detector
detector = prmflow.system.detector;

[fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice*Nfocal);

% equal angles
delta_fan = atan(mod(detector.mid_U, 1)*detector.hx_ISO/detector.SID) + ...
            atan((1-mod(detector.mid_U, 1))*detector.hx_ISO/detector.SID);
equalfan = ((1:Npixel)' - detector.mid_U).*delta_fan;
equalfan = repmat(equalfan + focalangle, Nslice, 1);
equalfan = reshape(equalfan, Npixel, Nslice*Nfocal);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nfocal, []);

% interp1
for ii = 1:Nslice*Nfocal
    [index_ii, alpha_ii] = interpprepare(fanangles(:, ii), equalfan(:, ii), 'extrap');
    dataflow.rawdata(:, ii, :) = dataflow.rawdata(index_ii(1,:), ii, :).*alpha_ii(1, :) + ...
                                 dataflow.rawdata(index_ii(2,:), ii, :).*alpha_ii(2, :);
end

% to return
prmflow.recon.fanangles = reshape(equalfan, Npixel*Nslice, Nfocal);
prmflow.recon.focalangle = focalangle;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end