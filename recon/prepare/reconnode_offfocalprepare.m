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
% do not use it in deployments
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
if isfield(prmflow.system, 'slicezebra')
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
% to ignore the difference of the focalspots in off-focal correction
focalposition = prmflow.system.focalposition(1, :);

% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice);

% % to ignore the error of detector position
% fanangles = mean(fanangles, 2);
% % We may save the fanangles in corr table

% slicemerge (in correction)
if ~isfield(offcorr, 'offslicemergescale')
    offcorr.offslicemergescale = 1;
end
offslicemergescale = offcorr.offslicemergescale;
Nslicemerged = Nslice/offcorr.offslicemergescale;

% merge the fanangles
if offslicemergescale > 1
    % slice merge the fanangles
    if ~slicezebra
        fanangles = squeeze(mean(reshape(fanangles, Npixel, offslicemergescale, Nslice/offslicemergescale), 2));
    else
        fanangles = mean(reshape(fanangles, Npixel, 2, offslicemergescale, Nslice/offslicemergescale/2), 3);
        fanangles = reshape(fanangles, Npixel, Nslice/offslicemergescale);
    end
end
% offcorr.Nslicemerged = Nslice/offslicemergescale;
% Note: the offslicemergescale is not the mergescale which is a general value in corr table for collimatorexposure. And the
% offslicemergescale is only for off-focal correction.

% viewsparse
if ~isfield(offcorr, 'viewsparse')
    offcorr.viewsparse = 1;
end
viewsparse = offcorr.viewsparse;

delta_view = pi*2 / Nviewprot * Nfocal;

% z-cross
if ~isfield(offcorr, 'crossrate')
    offcorr.crossrate = 0;
    % 0: the mean of all the slices
end
crossrate = offcorr.crossrate;

% z-cross matrix
crsMatrix = offfocalzcrossmatrix(Nslice, crossrate, slicezebra);
% the z-cross matrix works on not merged data

% Noffsample
if ~isfield(offcorr, 'Noffsample')
    offcorr.Noffsample = max(2^ceil(log2(Npixel)), 64);  % =1024
end
Noffsample = double(offcorr.Noffsample);    % int

% concyclic or homocentric
if isfield(detector, 'concyclic')
    isconcyclic = detector.concyclic;
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
    dt = offcorr.deltat/SID;
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
    airrate = squeeze(mean(reshape(airrate, Npixel, offcorr.offslicemergescale, Nslice/offcorr.offslicemergescale), 2));
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
% rawinterp2t = interp1(Dphi, (1:Npixel)', t_resp, 'linear', 'extrap');
% tinterp2raw is used to interp the off-focal fix to raw data space
tinterp2raw = zeros(Npixel, Nslicemerged, 'single');
% tinterp2raw = interp1(t_resp, (1:Noffsample)', Dphi, 'linear', 'extrap');
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
if ~isfield(offcorr, 'minintensity')
    if isfield(prmflow.raw, 'maxprojection')
        % maxprojection is about the photon-starvation correction (we don't have this value yet)
        offcorr.minintensity = 2^(-prmflow.raw.maxprojection);
        % plz ask a requirement to the photon-starvation calibration which need a standard condition
    else
        offcorr.minintensity = single(2^-32);
    end
end
% mincut of correction
if ~isfield(offcorr, 'fixcut')
    offcorr.fixcut = single(0.2);
end

% put these variables to corrtable
offcorr.crsMatrix = crsMatrix;
offcorr.slicezebra = slicezebra;
offcorr.delta_view = delta_view;
offcorr.offkernel = offkernel;
offcorr.Dphiscale = Dphiscale;
offcorr.Dphiscale_odd = Dphiscale_odd;
offcorr.Dphi = Dphi;
offcorr.rawinterp2t = rawinterp2t;
offcorr.rawinterp2phi = rawinterp2phi;
offcorr.Dphi = Dphi;
offcorr.tinterp2raw = tinterp2raw;
offcorr.tinterp2phi = tinterp2phi;
prmflow.corrtable.(status.nodename) = offcorr;
% we should update the off-focal calibration table to include those things

% These are only used while pipeline off:
prmflow.pipe.(nodename).offviewextra = offviewextra;
prmflow.pipe.(nodename).offstartview = offstartview;
prmflow.pipe.(nodename).offendview = offendview;
prmflow.pipe.(nodename).Nviewoff = Nviewoff;
% record Nshot and realated parameters in nodeprm (also, only used while pipeline off)
prmflow.pipe.(nodename).Nshot = prmflow.raw.Nshot;
prmflow.pipe.(nodename).viewpershot = prmflow.raw.viewpershot(...
    prmflow.raw.startshot : prmflow.raw.startshot + prmflow.raw.Nshot - 1);
% I worry some following nodes could change the Nshot

% pipe line
if pipeline_onoff
    % pipeline console paramters
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
    % viewcommon
    prmflow.pipe.(nodename).pipeline.viewcommon = viewsparse * Nfocal;

    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end

    % inner buffer
    dataflow.buffer.(nodename) = struct();
    % inner offspacepool
    dataflow.buffer.(nodename).offspacepool = status.defaultpool;
    % offspacepool.data and offspacepool.datafields
    dataflow.buffer.(nodename).offspacepool.data = struct();
    dataflow.buffer.(nodename).offspacepool.data.rawdata = single([]);
    dataflow.buffer.(nodename).offspacepool.datafields = {'rawdata'};


    % inner pipeline console
    prmflow.pipe.(nodename).raw2off.pipeline = defaultpipeprm();
    prmflow.pipe.(nodename).off2raw.pipeline = defaultpipeprm();
%     dataflow.buffer.(nodename).pipeline.raw2off = defaultpipeprm();
%     dataflow.buffer.(nodename).pipeline.off2raw = defaultpipeprm();  

    % inner raw2off 
    prmflow.pipe.(nodename).raw2off.pipeline.kernellevel = 1;
    prmflow.pipe.(nodename).raw2off.pipeline.viewrely = [offviewrely(2), offviewrely(1)].*viewsparse;
    prmflow.pipe.(nodename).raw2off.pipeline.viewrely_out = [offviewrely(2), offviewrely(1)];
    prmflow.pipe.(nodename).raw2off.pipeline.viewextra = offviewextra;
    prmflow.pipe.(nodename).raw2off.pipeline.viewrescale = [1 viewsparse];
    prmflow.pipe.(nodename).raw2off.pipeline.viewexpand = 0;
    prmflow.pipe.(nodename).raw2off.pipeline.viewcommon = viewsparse * Nfocal;
    prmflow.pipe.(nodename).raw2off.pipeline.inputminlimit = 1;
    prmflow.pipe.(nodename).raw2off.pipeline.inputmaxlimit = inf;
    prmflow.pipe.(nodename).raw2off.pipeline.iscarried = false;
    % inner off2raw
    prmflow.pipe.(nodename).off2raw.pipeline.kernellevel = 1;
    prmflow.pipe.(nodename).off2raw.pipeline.viewrely = offviewrely;
    prmflow.pipe.(nodename).off2raw.pipeline.viewrely_out = offviewrely.*viewsparse;
    prmflow.pipe.(nodename).off2raw.pipeline.viewextra = -offviewextra.*viewsparse;
    prmflow.pipe.(nodename).off2raw.pipeline.viewrescale = [viewsparse 1];
    prmflow.pipe.(nodename).off2raw.pipeline.viewexpand = 0;
    prmflow.pipe.(nodename).off2raw.pipeline.viewcommon = Nfocal;
    prmflow.pipe.(nodename).off2raw.pipeline.iscarried = false;
    prmflow.pipe.(nodename).off2raw.pipeline.inputmaxlimit = inf;
    
    % switch scan 
    switch lower(prmflow.protocol.scan)
        case 'axial'
            % the offspacepool is in size Nviewprot / viewsparse
            dataflow.buffer.(nodename).offspacepool.poolsize = prmflow.raw.Nviewprot / viewsparse;
            dataflow.buffer.(nodename).offspacepool.circulatemode = true;
            % inner raw2off 
            prmflow.pipe.(nodename).raw2off.pipeline.nodetype = 'A-A.1.G';
            prmflow.pipe.(nodename).raw2off.pipeline.relystrategy = 1;
            % inner off2raw
            prmflow.pipe.(nodename).off2raw.pipeline.nodetype = 'A-A.1.S';
            prmflow.pipe.(nodename).off2raw.pipeline.relystrategy = 2;
            prmflow.pipe.(nodename).off2raw.pipeline.inputminlimit = 1;
            % nextcirculte
            prmflow.pipe.(nodename).pipeline.nextcirculte = true;
        case {'helical', 'halfaxial', 'static', 'conveyor'}
            if ~isavail(dataflow.buffer.(nodename).offspacepool.poolsize)
                setpoolsize = prmflow.pipe.(nodename).pipeline.currpoolsize;
                % the offspacepool shall big enough
                if isavail(setpoolsize) % most not, unless user set it
                    dataflow.buffer.(nodename).offspacepool.poolsize = ...
                        ceil(setpoolsize / viewsparse) + sum(prmflow.pipe.(nodename).pipeline.viewrely) / viewsparse;
                else
                    dataflow.buffer.(nodename).offspacepool.poolsize = ...
                        ceil(prmflow.system.defaultrawpoolsize / viewsparse) + ...
                        sum(prmflow.pipe.(nodename).pipeline.viewrely) / viewsparse;
                end
                dataflow.buffer.(nodename).offspacepool.circulatemode = false;
            end
            % inner raw2off 
            prmflow.pipe.(nodename).raw2off.pipeline.nodetype = 'H-H.1.G';
            prmflow.pipe.(nodename).raw2off.pipeline.relystrategy = 1;
            % inner off2raw
            prmflow.pipe.(nodename).off2raw.pipeline.nodetype = 'H-H.1.S';
            prmflow.pipe.(nodename).off2raw.pipeline.relystrategy = 2;
            prmflow.pipe.(nodename).off2raw.pipeline.inputminlimit = -inf;  % force to ignore the minlimit
        otherwise
            0;
    end
end

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    % put corrtable to GPU
    prmflow.corrtable.(nodename) = putinGPU(prmflow.corrtable.(nodename));
%     % put parameters to GPU
%     GPUfields = {'crossrate', 'crsMatrix', 'offkernel', 'Dphiscale', 'Dphiscale_odd', 'Dphi', ...
%         'rawinterp2t', 'rawinterp2phi', 'tinterp2raw', 'tinterp2phi', 'minintensity'};
%     prmflow.correction.offfocal = putfieldsinGPU(prmflow.correction.offfocal, GPUfields);
end

% ini buffer..rawdata
if pipeline_onoff
    if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
        dataflow.buffer.(nodename).offspacepool.data.rawdata = ...
            zeros(Noffsample*Nslicemerged, dataflow.buffer.(nodename).offspacepool.poolsize, 'single', 'gpuArray');
    else
        dataflow.buffer.(nodename).offspacepool.data.rawdata = ...
            zeros(Noffsample*Nslicemerged, dataflow.buffer.(nodename).offspacepool.poolsize, 'single');
    end
end


end