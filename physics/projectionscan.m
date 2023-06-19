function Dataflow = projectionscan(SYS, method, echo_onoff)
% the projection simulation

if nargin<2 || isempty(method)
    method = 'default';
end
if nargin<3
    echo_onoff = true;
end

% GPU?
GPUonoff = SYS.simulation.GPUonoff>0 && ~isempty(SYS.simulation.GPUinfo);

% system components
source = SYS.source;
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
if isfield(SYS, 'phantom')
    phantom = SYS.phantom;
else
    phantom = [];
end

% parameters of the system
% focalposition = source.focalposition;
% Nfocal = source.focalnumber;
Nfocalpos = size(source.focalposition, 1);
Npixel = double(detector.Npixel);
Nslice = double(detector.Nslice);
Np = Npixel * Nslice;
Nw = source.Wnumber;

% cut pixels?
if isfield(detector, 'pixelrange') && ~isempty(detector.pixelrange)
    pixelrange = double(detector.pixelrange);
    Np = double(detector.Nprange)*Nslice;
else
    pixelrange = [];
end
if ~isempty(pixelrange)
    Npr = size(pixelrange, 2);
    if Npr ~= Nfocalpos
        pixelrange = repelem(pixelrange, 1, Nfocalpos/Npr);
    end
end

% prepare the samplekeV, viewangle and couch
[samplekeV, viewangle, couch, shotindex, gantrytilt] = scanprepare(SYS);
NkeVsample = length(samplekeV(:));
Nview = length(viewangle(:));
Nviewpf = Nview/Nfocalpos;

% spectrums normalize
sourcespect = SYS.source.spectrum;
for ii = 1:Nw
    sourcespect{ii} = sourcespect{ii}./sum(sourcespect{ii}.*samplekeV);
end
% detector response
detspect = cell(1, Nw);
respnorm = sum(samplekeV)/mean(detector.response * samplekeV(:), 1);
for ii = 1:Nw
    detspect{ii} = detector.response.*sourcespect{ii}.*respnorm;
    % if only one detspect curve, do nothing
    if size(detspect{ii}, 1) ~= 1
        if isempty(pixelrange)
            detspect{ii} = repmat(detspect{ii}, Np*Nfocalpos, 1);
        else
            detspect_tmp = detspect{ii};
            detspect{ii} = zeros(Np*Nfocalpos, NkeVsample, class(detspect_tmp));
            for ifocal = 1:Nfocalpos
                index_ifoc = (1:Np) + (ifocal-1)*Np;
                index_det = mod(pixelrange(1,ifocal)-1+(0:Nprange-1)', Npixel)+1 + (0:Nslice-1).*Npixel;
                detspect{ii}(index_ifoc, :) = detspect_tmp(index_det(:), :);
            end
        end
    end
end
% detector pixelarea
detpixelarea = detector.pixelarea(:);
% if size(detector.pixelarea(:), 1) == 1
%     detpixelarea = detector.pixelarea;
% else
%     detpixelarea = repmat(detector.pixelarea(:), 1, Nfocalpos, 1);
% end

% detector norminal vector
detnormv = detector.normvector;

% detector resample
if any(SYS.simulation.detectsample>1)
    [detposition, detweight, detwgrid] = detectorresample(detector, SYS.simulation.detectsample);
else
    detposition = detector.position;
    detweight = 1;
    detwgrid = [];
end
Ndetresmp = size(detweight, 1);

% focal resample
if any(SYS.simulation.focalsample>1)
    [focalposition, focalweight] = focalresample(source.focalposition, source.focalsize, SYS.simulation.focalsample);
else
    focalposition = source.focalposition;
    focalweight = 1;
end
Nfocalresmp = size(focalweight, 1);

% off-focal
if SYS.simulation.offfocal && isfield(SYS.source, 'offfocal')
    [offfocalpos, offfocalweight] = offfocalresample(source.focalposition, SYS.source.offfocal);
    Nofffocalsmp = size(offfocalweight, 1);
else
    offfocalpos = [];
    offfocalweight = [];
    Nofffocalsmp = 0;
end
focalposition = [focalposition offfocalpos];

% det-focal resample grid
[focalindex, detindex] = meshgrid(1:Nfocalresmp+Nofffocalsmp, 1:Ndetresmp);
weight = detweight*[focalweight; offfocalweight]';
Nresample = size(weight(:),1);

% ASG
ASGdata = prepareasgdata(detector, samplekeV);

% ini P, Pair  & Eeff
P = cell(1, Nw);
Pair = cell(1, Nw);
Eeff = cell(1, Nw);
switch lower(method)
    case {'default', 1, 'photoncount', 2}
        P(:) = {zeros(Np, Nview)};
        Pair(:) = {zeros(Np*Nfocalpos, 1)};
        Eeff(:) = {zeros(Np, Nview)};
    case {'energyvector', 3}
        P(:) = {zeros(Np*Nview, NkeVsample)};
        Pair(:) = {zeros(Np*Nfocalpos, NkeVsample)};
        % No Eeff
    otherwise
        % error
        error(['Unknown projection method: ' method]);
end


for ismp = 1:Nresample
    focalposition_ismp = focalposition(:, (1:3)+(focalindex(ismp)-1)*3);
    detposition_ismp = detposition(:, (1:3)+(detindex(ismp)-1)*3);
    if ~isempty(ASGdata) && ~isempty(detwgrid)
        % gap resample
        ASGdata = gapresample(ASGdata, detwgrid(detindex(ismp), :));
    end
    % I know the detnormv of a detector pxiel is not change after any resampling
    [P_ismp, Pair_ismp, Eeff_ismp] = projectionscan2(focalposition_ismp, detposition_ismp, Npixel, pixelrange, detnormv, ...
        bowtie, filter, samplekeV, detspect, detpixelarea, viewangle, couch, gantrytilt, phantom, method, ASGdata, ...
        echo_onoff, GPUonoff);
    for iw = 1:Nw
        P{iw} = P{iw} + P_ismp{iw}.*weight(ismp);
        Pair{iw} = Pair{iw} + Pair_ismp{iw}.*weight(ismp);
        Eeff{iw} = Eeff{iw} + Eeff_ismp{iw}.*weight(ismp);
    end
end
for iw = 1:Nw
    Eeff{iw} = Eeff{iw}./sum(weight(:));
end

% remove nan (could due to negative thick filter)
for ii = 1:Nw
    P{ii}(isnan(P{ii})) = 0;
end

% return
Dataflow = struct();
Dataflow.samplekeV = samplekeV;
Dataflow.viewangle = viewangle;
Dataflow.couch = couch;
Dataflow.shotindex = shotindex;
Dataflow.P = P;
Dataflow.Eeff = Eeff;
Dataflow.Pair = Pair;

end


function [v, w, wgrid] = detectorresample(detector, Nresmp)
% detector resample

Nx = Nresmp(1);
if size(Nresmp(:),1)>1
    Nz = Nresmp(2);
else
    Nz = 1;
end

% sample nodes and weight
[px, wx] = lgwt(Nx, -0.5, 0.5);
[pz, wz] = lgwt(Nz, -0.5, 0.5);
w = wx*wz';
w = w(:);
Np = size(w, 1);
[PZ, PX] = meshgrid(pz,px);
% wgrid
gwx = [0; cumsum(wx)];
gwz = [0; cumsum(wz)];
wgrid =[repmat([gwx(1:end-1) 1-gwx(2:end)], Nz, 1)  repelem([gwz(1:end-1) 1-gwz(2:end)], Nx, 1)];

% x
if isfield(detector, 'nVedgex')
    Vx = detector.nVedgex;
else
    Vx = zeros(size(detector.normvector));
    Vx(:, 1) = detector.normvector(:, 2);
    Vx(:, 2) = -detector.normvector(:, 1);
    Vx = normr(Vx);
end
% z
if isfield(detector, 'nVedgez')
    Vz = detector.nVedgez;
else
    Vz = cross(Vx, detector.normvector);
end
% length scale
Vx = Vx.*detector.edgelength(:, 1);
Vz = Vz.*detector.edgelength(:, 2);

% v
v = detector.position(:) + (Vx(:)*PX(:)' + Vz(:)*PZ(:)');
v = reshape(v, [], 3*Np);

end


function [v, w] = focalresample(focalpos, focalsize, Nresmp)
% detector resample

Nx = Nresmp(1);
if size(Nresmp(:),1)>1
    Nz = Nresmp(2);
else
    Nz = 1;
end
Nfocal = size(focalpos, 1);

% sample nodes and weight
[px, wx] = lgwt(Nx, -0.5, 0.5);
[pz, wz] = lgwt(Nz, -0.5, 0.5);
w = wx*wz';
w = w(:);
Np = size(w, 1);
[PZ, PX] = meshgrid(pz,px);
% Vx Vz
Vx = repmat([focalsize(1) 0 0], Nfocal, 1);
Vz = repmat([0 0 focalsize(2)], Nfocal, 1);

% v
v = focalpos(:) + (Vx(:)*PX(:)' + Vz(:)*PZ(:)');
v = reshape(v, [], 3*Np);

end


function [offfocalpos, offfocalweight] = offfocalresample(focalposition, offfocal)
% off-focal resample

Nfocal = size(focalposition, 1);

if isfield(offfocal, 'samplepos')
    offfocalpos = focalposition(:) + repelem(reshape(offfocal.samplepos, 3, []), Nfocal, 1);
    offfocalpos = reshape(offfocalpos, Nfocal, []);
    offfocalweight = offfocal.weight(:);
elseif isfield(offfocal, 'width')
    % simplified defination of off-focal (hard code)
    Nsf = 10;
    X = (-Nsf+1/2:Nsf-1/2).*(offfocal.width/Nsf);
    offfocalpos = focalposition(:) + [repelem(X, Nfocal, 1); zeros(Nfocal*2, Nsf*2)];
    offfocalpos = reshape(offfocalpos, Nfocal, []);
    offfocalweight = ones(Nsf*2, 1).*(offfocal.intensity/Nsf/2);
else
    offfocalpos = [];
    offfocalweight = [];
end

end


function ASGdata = gapresample(ASGdata, detwgrid)
% resample of ASG gap
edgegapX = ASGdata.edgegap(:, 1:2) + ASGdata.edgelength(:, 1)*detwgrid(1:2);
edgegapZ = ASGdata.edgegap(:, 3:4) + ASGdata.edgelength(:, 2)*detwgrid(3:4);

ASGdata.gap_p = [(edgegapX(:, 1) + edgegapX(:, 2))./2  (edgegapZ(:, 1) + edgegapZ(:, 2))./2];
ASGdata.gap_n = [(edgegapX(:, 1) - edgegapX(:, 2))./2  (edgegapZ(:, 1) - edgegapZ(:, 2))./2];
end