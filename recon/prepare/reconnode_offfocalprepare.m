function [dataflow, prmflow, status] = reconnode_offfocalprepare(dataflow, prmflow, status)
% recon node, off-focal correction prepare
% [dataflow, prmflow, status] = reconnode_offfocalprepare(dataflow, prmflow, status);

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

% parameters of this node in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    offcorr = prmflow.corrtable.(status.nodename);
else
    error('Did not load corrtable for %s!', status.nodename);
end
offcorrversion = str2double([num2str(offcorr.ID(3)) '.' num2str(offcorr.ID(4))]);

% off-focal kernel parameter (used in cali, in manually setting off-focalkernel and/or simulation off-focal)
if isfield(nodeprm, 'offfocalkernel')
    offfocalkernel = nodeprm.offfocalkernel;
elseif isfield(prmflow.system, 'offfocalkernel')
    offfocalkernel = prmflow.system.offfocalkernel;
else
    offfocalkernel = [];
end
if ~isempty(offfocalkernel)
    % load offfocalkernel file
    if ischar(offfocalkernel)
        offfocalkernel = readcfgfile(offfocalkernel);
    end
    offcorr_fromkernel = offfocalloadkernel(offfocalkernel, prmflow.protocol);
    % merge the offcorr_fromkernel and/or caliprm to offcorr
    nodeprm = structmerge(nodeprm, offcorr_fromkernel, 1, 0); 
end
% merge the parapmeters to offcorr
offcorr = structmerge(nodeprm, offcorr, 1, 0);

% slice cross designment
if isfield(offcorr, 'slicezebra')
    slicezebra = offcorr.slicezebra;
elseif isfield(prmflow.system, 'slicezebra')
    slicezebra = prmflow.system.slicezebra;
else
    slicezebra = false;
end

% parameters to use in prmflow
Nshot = prmflow.raw.Nshot;
Nview = prmflow.raw.Nview;
Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
Nviewprot = prmflow.raw.Nviewprot;
scantype = prmflow.raw.scan;

% detector
detector = prmflow.system.detector;
SID = detector.SID;
SDD = detector.SDD;
% focalspot = prmflow.raw.focalspot;
% ignore the difference of the focalspots in off-focal correction
focalposition = prmflow.system.focalposition(1, :);

% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice);

% slicemerge
if isfield(nodeprm, 'slicemerge')
    slicemerge = nodeprm.slicemerge;
    if ~slicezebra
        fanangles = squeeze(mean(reshape(fanangles, Npixel, slicemerge, Nslice/slicemerge), 2));
    else
        fanangles = mean(reshape(fanangles, Npixel, 2, slicemerge, Nslice/slicemerge/2), 3);
        fanangles = reshape(fanangles, Npixel, Nslice/slicemerge);
    end
else
    slicemerge = 1;
end
Nslicemerge = Nslice/slicemerge;

% viewsparse
if isfield(offcorr, 'viewsparse')
    viewsparse = offcorr.viewsparse;
else
    viewsparse = 1;
end
delta_view = pi*2 / Nviewprot;

% z-cross
if isfield(offcorr, 'crossrate')
    if isfinite(offcorr.crossrate)
        crossrate = offcorr.crossrate;
    else
        warning('Illeagal value of ''crossrate'' in off-focal calibration table! Replaced by 0.');
        crossrate = 0;
    end
else
    crossrate = 0;
    % 0 is the mean of all the slices
end
% z-cross matrix
if ~slicezebra
    crsMatrix = offfocalzcrossmatrix(Nslice, crossrate);
else
    crsMatrix = offfocalzcrossmatrix(Nslice/2, crossrate);
end

% Noffsample
if isfield(offcorr, 'offsample')
    Noffsample = nodeprm.offsample;
    % To set the Noffsample in 1024 or 512
else
    Noffsample = max(2^ceil(log2(Npixel)), 64);  % =1024
end

% off-focal tau-measure
phi = fanangles-pi/2;
alpha = acos(SID/SDD);
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);

% extra view
extraview = [floor(double(min(phi(:)) - max(phi_off(:)))/delta_view/viewsparse); ...
    ceil(double(max(phi(:)) - min(phi_off(:)))/delta_view/viewsparse)];

% start-end view of off-focal
switch lower(scantype)
    case {'axial', 'static'}
        offstartview = 1;
        offendview = Nviewprot/viewsparse;   % must be integer
        Nviewoff = Nviewprot/viewsparse;
    case {'helical', 'halfaxial'}
        if isfield(nodeprm, 'startend_flag')
            startend_flag = nodeprm.startend_flag;
        else
            startend_flag = 0;
        end
        Nviewsparse = ceil(Nview/viewsparse);
        switch startend_flag
            case 1
                % middle
                offstartview = 1 + floor(double(-max(phi_off(:)))/delta_view/viewsparse);
                offendview = Nviewsparse + ceil(double(-min(phi_off(:)))/delta_view/viewsparse);
            case 2
                % full correction
                offstartview = 1 + extraview(1);
                offendview = Nviewsparse + extraview(2);
            otherwise
                % simplified
                offstartview = 1;
                offendview = Nviewsparse;
        end
        Nviewoff = offendview - offstartview + 1;
    otherwise
        % what?
        offstartview = 0;
        offendview = -1;
        Nviewoff = 0;
end

% Dphi and off-focal sampling
Dphi = phi - phi_off;
% Dphi_mean = mean(Dphi, 2);

% remeasurement, t0
if offcorrversion < 2.0
    % old version 1.x
    % the Riemannian measure is tan/atan, hard coded.
    minDphi = tan(min(Dphi(:)));
    maxDphi = tan(max(Dphi(:)));
else
    minDphi = min(Dphi(:));
    maxDphi = max(Dphi(:));
end
% dt and t0, dt = deltat is configurable.
if isfield(offcorr, 'deltat')
    dt = deltat/SID;
    t0 = (0 : Noffsample-1)' .* dt + minDphi;
    if max(t0) < maxDphi
        warning('The parameter offsample is too small!');
    end
else
    dt = (maxDphi - minDphi)/(Noffsample - 1);
    t0 = linspace(minDphi, maxDphi, Noffsample)';
end
% remeasurement, t_resp
if offcorrversion < 2.0
    % tan/atan measure
    t_resp = atan(t0);
    % the Dphiscale was airrate
    airrate = 2.^(-reshape(offcorr.airrate, Npixel, Nslice).*offcorr.ratescale(1));
    airrate = squeeze(mean(reshape(airrate, Npixel, slicemerge, Nslice/slicemerge), 2));
    Dphiscale = max(airrate, 1)./airrate;
    Dphiscale_odd = 0;
else
    % numerical measure
    t_resp = interp1(offcorr.t0, offcorr.tresp, t0, 'linear', 'extrap');
    Dphiscale = interp1(offcorr.t0, offcorr.tscale, Dphi, 'linear', 'extrap');
    Dphiscale_odd = interp1(offcorr.t0, offcorr.tscaleodd, Dphi, 'linear', 'extrap');
end
% rawinterp2t is used to interp the raw data to off-focal measurment space along the channel direction
rawinterp2t = zeros(Noffsample, Nslicemerge, 'single');
% tinterp2raw is used to interp the off-focal fix to raw data space
tinterp2raw = zeros(Npixel, Nslicemerge, 'single');
for ii = 1:Nslicemerge
    rawinterp2t(:, ii) = interp1(Dphi(:, ii), (1:Npixel)', t_resp, 'linear', 'extrap');
    tinterp2raw(:, ii) = interp1(t_resp, (1:Noffsample)', Dphi(:, ii), 'linear', 'extrap');
end
rawinterp2phi = t_resp./delta_view;
tinterp2phi = Dphi./delta_view/viewsparse;

% off-focal kernel
if offcorrversion < 2.0
    % old version 1.x
    % sinc kernel
    offwidth_nrm = offcorr.offwidth/SID/(maxDphi - minDphi);
    offkernel = offfocalsinckernel(offcorr.offintensity, offwidth_nrm, offcorr.offedge, Noffsample);
    % normalization by setting:
    offkernel = offkernel - offkernel(1);
else
    % numerical off-focal kernel
    % resample
    dx = SID*dt;
    Ncorrsamp = (length(offcorr.x)-1)/2;
    Nkernelsamp = ceil((dx/offcorr.dx) * Ncorrsamp);
    xksamp = (-Nkernelsamp : Nkernelsamp)' .* dx;
    curve_corr = (offcorr.curve + offcorr.curveodd.*offcorr.Iodd.*1i).*offcorr.Iscale;
    offcurve = spectrumresample(offcorr.x, curve_corr, xksamp);
    g2s = zeros(Noffsample, 1, 'single');
    g2s(1:Nkernelsamp+1) = offcurve(Nkernelsamp+1 : end);
    g2s(end-Nkernelsamp+1:end) = offcurve(1:Nkernelsamp);
    offkernel = 1./(1-fft(g2s))-1;
    % normalization
    offkernel = offkernel - offkernel(1);
end

% minimum intensity
if isfield(offcorr, 'minintensity')
    minintensity = offcorr.minintensity;
elseif isfield(prmflow.raw, 'maxprojection')
    % maxprojection is about the photon-starvation correction
    minintensity = 2^(-prmflow.raw.maxprojection);
else
    minintensity = 2^-32;
end

% save to prmflow, in .raw.offfocal
prmflow.raw.offfocal.crossrate = crossrate;
prmflow.raw.offfocal.crsMatrix = crsMatrix;
prmflow.raw.offfocal.slicezebra = slicezebra;
prmflow.raw.offfocal.slicemerge = slicemerge;
prmflow.raw.offfocal.viewsparse = viewsparse;
prmflow.raw.offfocal.delta_view = delta_view;
prmflow.raw.offfocal.offstartview = offstartview;
prmflow.raw.offfocal.offendview = offendview;
prmflow.raw.offfocal.Nviewoff = Nviewoff;
prmflow.raw.offfocal.extraview = extraview;
prmflow.raw.offfocal.Noffsample = Noffsample;
prmflow.raw.offfocal.offkernel = offkernel;
prmflow.raw.offfocal.Dphiscale = Dphiscale;
prmflow.raw.offfocal.Dphiscale_odd = Dphiscale_odd;
prmflow.raw.offfocal.Dphi = Dphi;
prmflow.raw.offfocal.rawinterp2t = rawinterp2t;
prmflow.raw.offfocal.rawinterp2phi = rawinterp2phi;
prmflow.raw.offfocal.tinterp2raw = tinterp2raw;
prmflow.raw.offfocal.tinterp2phi = tinterp2phi;
prmflow.raw.offfocal.minintensity = minintensity;

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpooldata;
    dataflow.buffer.(nodename) = struct();
    % output pool
    dataflow.buffer.(nodename).outpool = struct();
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = 1;
    dataflow.buffer.(nodename).AvailViewindex = 0;
    dataflow.buffer.(nodename).AvailPoint = 0;
    dataflow.buffer.(nodename).outpoolsize = Inf;
    % inner buffer
    dataflow.buffer.(nodename).offfocalspace = struct();
    dataflow.buffer.(nodename).offReadPoint = 1;
    dataflow.buffer.(nodename).offWritePoint = 1;
    dataflow.buffer.(nodename).offReadViewindex = offstartview;
    dataflow.buffer.(nodename).offpoolsize = Inf;
end

end