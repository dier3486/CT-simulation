function [dataflow, prmflow, status] = reconnode_detectorcali(dataflow, prmflow, status)
% detector calibration
% [dataflow, prmflow, status] = reconnode_detectorcali(dataflow, prmflow, status)

% parameters and data to use in prmflow
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nviewprot = prmflow.recon.Nviewprot;
detector = prmflow.system.detector;
focalpos = prmflow.system.focalposition(prmflow.recon.focalspot, :);

% parameters to from pipe
caliprm = prmflow.pipe.(status.nodename);

% Npixelpermod
if isfield(caliprm, 'Npixelpermod')
    Npixelpermod = caliprm.Npixelpermod;
elseif isfield(detector, 'Npixelpermod')
    Npixelpermod = detector.Npixelpermod;
else
    % default
    Npixelpermod = 16;
end

% pinfit ini &  fitting coefficients
if isfield(caliprm, 'pinfit_ini')
    pinfit_ini = caliprm.pinfit_ini;
else
%                [r,       phi,     zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, isooff2, isophase2, zshift]
    pinfit_ini = [0.4      pi/2     0       0       0         0     0       0       0       0];
end
if isfield(caliprm, 'fitselect')
    fitselect = caliprm.fitselect;
else
    fitselect =  [1        1        1       1       0         1     1       1       1       1];
end
% I know the pinfit(2) is the angle of the pin located, set the pinfit_ini(2) to p/2 when the pin is below the ISO.
% I know the pinfit(5) is the 'rotscale' to used when the pin can not matched in a whole rotation, e.g. a time based scan.
% Please set the fitselect(5) to false to skip that function in angle based scans.

% alphaL
if isfield(caliprm, 'alphaL')
    alphaL = caliprm.alphaL;
else
    alphaL = [1.5 0.1];
end
% fanfix
[fanfix, pinfit] = pindetcalibration(dataflow.rawdata, dataflow.rawhead.viewangle, detector, focalpos, ...
                   Npixel, Nslice, Npixelpermod, Nviewprot, alphaL, pinfit_ini, fitselect);

% apply on XYZ 
midfanfix = pinfit(6);  % I know the pinfit(6) is the fixing of mid channel
Xfix = (detector.position(:, 1) - focalpos(1)).*cos(fanfix(:) - midfanfix) - ...
       (detector.position(:, 2) - focalpos(2)).*sin(fanfix(:) - midfanfix) + focalpos(1);
Yfix = (detector.position(:, 1) - focalpos(1)).*sin(fanfix(:) - midfanfix) + ...
       (detector.position(:, 2) - focalpos(2)).*cos(fanfix(:) - midfanfix) + focalpos(2);
    
detectorcorr = detector;
detectorcorr.position(:, [1 2]) = [Xfix, Yfix];

% to return
dataflow.detectorcorr = detectorcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
