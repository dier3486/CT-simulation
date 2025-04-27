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
Nfocal = prmflow.raw.Nfocal;
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
Nslicemerged = Nslice/slicemerge;

% viewsparse
if isfield(offcorr, 'viewsparse')
    viewsparse = offcorr.viewsparse;
else
    viewsparse = 1;
end
delta_view = pi*2 / Nviewprot * Nfocal;

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
crsMatrix = offfocalzcrossmatrix(Nslice, crossrate, slicezebra);

% Noffsample
if isfield(offcorr, 'Noffsample')
    Noffsample = single(offcorr.Noffsample);
    % To set the Noffsample in 1024 or 512
else
    Noffsample = max(2^ceil(log2(Npixel)), 64);  % =1024
end

% concyclic or homocentric
if isfield(offcorr, 'concyclic')
    isconcyclic = offcorr.concyclic;
elseif isfield(detector, 'concyclic')
    isconcyclic = detector.concyclic;
elseif isfield(nodeprm, 'concyclic')
    isconcyclic = nodeprm.concyclic;
else
    % default is homocentric (not concyclic-detector)
    isconcyclic = false;
end

% off-focal tau-measure
phi = fanangles-pi/2;
alpha = acos(SID/SDD);
if ~isconcyclic
    % homocentric
    phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);
else
    % concyclic
    phi_off = phi.*(SID/(SID-SDD));
end

% off-focal view rely
offviewrely = [ceil(double(max(phi_off(:)) - min(phi(:)))/delta_view/viewsparse), ...
    ceil(double(max(phi(:)) - min(phi_off(:)))/delta_view/viewsparse)] .* Nfocal;

% start-end view of off-focal

switch lower(scantype)
    case 'axial'
        offviewextra = [0 0];
        Nviewsparse = Nviewprot/viewsparse;
    case 'static'
        offviewextra = [0 0];
        Nviewsparse = floor(Nview/viewsparse);
    case {'helical', 'halfaxial'}
        if isfield(nodeprm, 'startend_flag')
            startend_flag = nodeprm.startend_flag;
        else
            startend_flag = 0;
        end
        Nviewsparse = floor(Nview/viewsparse);
        switch startend_flag
            case 1
                % middle
                offviewextra = ceil( double([max(phi_off(:)), -min(phi_off(:))])./delta_view./viewsparse ) .* Nfocal;
            case 2
                % full correction
                offviewextra = offviewrely;
            otherwise
                % simplified
                offviewextra = [0 0];
        end

    otherwise
        % what?
        offviewextra = [0 0];
        Nviewsparse = floor(Nview/viewsparse);
        % 
end
offstartview = 1 - offviewextra(1);
offendview = Nviewsparse + offviewextra(2);
Nviewoff = offendview - offstartview + 1;

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
rawinterp2t = zeros(Noffsample, Nslicemerged, 'single');
% tinterp2raw is used to interp the off-focal fix to raw data space
tinterp2raw = zeros(Npixel, Nslicemerged, 'single');
for ii = 1:Nslicemerged
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
    Ncorrsamp = (length(offcorr.curvesamp)-1)/2;
    Nkernelsamp = ceil((dx/offcorr.dcs) * Ncorrsamp);
    xksamp = (-Nkernelsamp : Nkernelsamp)' .* dx;
    curve_corr = (offcorr.curve + offcorr.curveodd.*offcorr.intensityodd.*1i).*offcorr.intensityscale;
    offcurve = spectrumresample(offcorr.curvesamp, curve_corr, xksamp);
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

% save to prmflow, in .correction.offfocal
prmflow.correction.offfocal.crossrate = crossrate;
prmflow.correction.offfocal.crsMatrix = crsMatrix;
prmflow.correction.offfocal.slicezebra = slicezebra;
prmflow.correction.offfocal.slicemerge = slicemerge;
prmflow.correction.offfocal.viewsparse = viewsparse;
prmflow.correction.offfocal.delta_view = delta_view;
prmflow.correction.offfocal.offstartview = offstartview;
prmflow.correction.offfocal.offendview = offendview;
prmflow.correction.offfocal.offviewextra = offviewextra;
prmflow.correction.offfocal.Nviewoff = Nviewoff;
prmflow.correction.offfocal.offviewrely = offviewrely;
prmflow.correction.offfocal.Noffsample = Noffsample;
prmflow.correction.offfocal.offkernel = offkernel;
prmflow.correction.offfocal.Dphiscale = Dphiscale;
prmflow.correction.offfocal.Dphiscale_odd = Dphiscale_odd;
prmflow.correction.offfocal.Dphi = Dphi;
prmflow.correction.offfocal.rawinterp2t = rawinterp2t;
prmflow.correction.offfocal.rawinterp2phi = rawinterp2phi;
prmflow.correction.offfocal.tinterp2raw = tinterp2raw;
prmflow.correction.offfocal.tinterp2phi = tinterp2phi;
prmflow.correction.offfocal.minintensity = minintensity;

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    % the off-focal correction is H-H.0.S or A.0.S 
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    if strcmpi(prmflow.protocol.scan, 'static')
        % but static scan is in type H.0.N
        prmflow.pipe.(nodename).pipeline.viewrely = [0 0];
        prmflow.pipe.(nodename).pipeline.relystrategy = 0;
    else
        viewrely = (offviewrely(1) + offviewrely(2)) * viewsparse;
        prmflow.pipe.(nodename).pipeline.viewrely = [viewrely viewrely];
        prmflow.pipe.(nodename).pipeline.relystrategy = 'stingy';
%         prmflow.pipe.(nodename).pipeline.inputminlimit = viewrely + 1;
    end
    prmflow.pipe.(nodename).pipeline.viewcommon = viewsparse * Nfocal;

    % inner buffer
    dataflow.buffer.(nodename) = struct();
    % inner offspacepool
    dataflow.buffer.(nodename).offspacepool = status.defaultpool;
    % inner pipeline console
    prmflow.correction.offfocal.pipeline = struct();
    prmflow.correction.offfocal.pipeline.raw2off = struct();
    prmflow.correction.offfocal.pipeline.off2raw = struct();

    if isfield(nodeprm, 'pipeline')
        % curr pool configured by user
        dataflow.pipepool.(nodename) = initialpool(dataflow.pipepool.(nodename), ...
            prmflow.pipe.(nodename).pipeline, 'public');
        % inner pool
        dataflow.buffer.(nodename).offspacepool = initialpool(dataflow.buffer.(nodename).offspacepool, ...
            prmflow.pipe.(nodename).pipeline, 'offspace');
    end
    % offspacepool.data and offspacepool.datafields
    dataflow.buffer.(nodename).offspacepool.data = struct();
    dataflow.buffer.(nodename).offspacepool.data.rawdata = single([]);
    dataflow.buffer.(nodename).offspacepool.datafields = {'rawdata'};
    % inner raw2off 
    prmflow.correction.offfocal.pipeline.raw2off.kernellevel = 1;
    prmflow.correction.offfocal.pipeline.raw2off.viewrely = [offviewrely(2), offviewrely(1)].*viewsparse;
    prmflow.correction.offfocal.pipeline.raw2off.viewrely_out = [offviewrely(2), offviewrely(1)];
    prmflow.correction.offfocal.pipeline.raw2off.viewextra = offviewextra;
    prmflow.correction.offfocal.pipeline.raw2off.viewrescale = [1 viewsparse];
    prmflow.correction.offfocal.pipeline.raw2off.viewexpand = 0;
    prmflow.correction.offfocal.pipeline.raw2off.viewcommon = viewsparse * Nfocal;
    prmflow.correction.offfocal.pipeline.raw2off.inputminlimit = 1;
    prmflow.correction.offfocal.pipeline.raw2off.inputmaxlimit = inf;
    prmflow.correction.offfocal.pipeline.raw2off.iscarried = false;
    % inner off2raw
    prmflow.correction.offfocal.pipeline.off2raw.kernellevel = 1;
    prmflow.correction.offfocal.pipeline.off2raw.viewrely = offviewrely;
    prmflow.correction.offfocal.pipeline.off2raw.viewrely_out = offviewrely.*viewsparse;
    prmflow.correction.offfocal.pipeline.off2raw.viewextra = -offviewextra.*viewsparse;
    prmflow.correction.offfocal.pipeline.off2raw.viewrescale = [viewsparse 1];
    prmflow.correction.offfocal.pipeline.off2raw.viewexpand = 0;
    prmflow.correction.offfocal.pipeline.off2raw.viewcommon = Nfocal;
    prmflow.correction.offfocal.pipeline.off2raw.iscarried = false;
    prmflow.correction.offfocal.pipeline.off2raw.inputmaxlimit = inf;
    
    % switch scan 
    switch lower(prmflow.protocol.scan)
        case 'axial'
            % the offspacepool is in size Nviewprot / viewsparse
            dataflow.buffer.(nodename).offspacepool.poolsize = prmflow.raw.Nviewprot / viewsparse;
            dataflow.buffer.(nodename).offspacepool.circulatemode = true;
            % inner raw2off 
            prmflow.correction.offfocal.pipeline.raw2off.nodetype = 'A-A.1.G';
            prmflow.correction.offfocal.pipeline.raw2off.relystrategy = 1;
            % inner off2raw
            prmflow.correction.offfocal.pipeline.off2raw.nodetype = 'A-A.1.S';
            prmflow.correction.offfocal.pipeline.off2raw.relystrategy = 2;
            prmflow.correction.offfocal.pipeline.off2raw.inputminlimit = 1;
            % nextcirculte
            prmflow.pipe.(nodename).pipeline.nextcirculte = true;
        case {'helical', 'halfaxial', 'static'}
            if ~isavail(dataflow.buffer.(nodename).offspacepool.poolsize)
                % the offspacepool shall big enough
                if isavail(dataflow.pipepool.(nodename).poolsize)
                    dataflow.buffer.(nodename).offspacepool.poolsize = ...
                        ceil(dataflow.pipepool.(nodename).poolsize / viewsparse) + ...
                        sum(prmflow.pipe.(nodename).pipeline.viewrely) / viewsparse;
                else
                    dataflow.buffer.(nodename).offspacepool.poolsize = ...
                        ceil(prmflow.system.defaultrawpoolsize / viewsparse) + ...
                        sum(prmflow.pipe.(nodename).pipeline.viewrely) / viewsparse;
                end
                dataflow.buffer.(nodename).offspacepool.circulatemode = false;
            end
            % inner raw2off 
            prmflow.correction.offfocal.pipeline.raw2off.nodetype = 'H-H.1.G';
            prmflow.correction.offfocal.pipeline.raw2off.relystrategy = 1;
            % inner off2raw
            prmflow.correction.offfocal.pipeline.off2raw.nodetype = 'H-H.1.S';
            prmflow.correction.offfocal.pipeline.off2raw.relystrategy = 2;
            prmflow.correction.offfocal.pipeline.off2raw.inputminlimit = -inf;  % force to ignore the minlimit
        otherwise
            0;
    end
    
    % ini buffer..rawdata
    dataflow.buffer.(nodename).offspacepool.data.rawdata = ...
        zeros(Noffsample*Nslicemerged, dataflow.buffer.(nodename).offspacepool.poolsize, 'single');
%     dataflow.buffer.(nodename).offspacepool.data.reading = ...
%         zeros(1, dataflow.buffer.(nodename).offspacepool.poolsize, 'single');
end


end