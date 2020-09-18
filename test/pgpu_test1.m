% function Dataflow = projectionscan(SYS, method, echo_onoff)
% the projection simulation

method = 'default';
echo_onoff = true;

% ini GPU Device
gpuD = gpuDevice;

tic;

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
focalposition = source.focalposition;
Nfocal = source.focalnumber;
Npixel = double(detector.Npixel);
Nslice = double(detector.Nslice);
Np = Npixel * Nslice;
Nw = source.Wnumber;
Mlimit = SYS.simulation.memorylimit;

% prepare the samplekeV, viewangle and couch
[samplekeV, viewangle, couch, shotindex, gantrytilt] = scanprepare(SYS);
Nsample = length(samplekeV(:));
Nview = length(viewangle(:));

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
    if size(detspect{ii}, 1) == 1
        detspect{ii} = repmat(detspect{ii}, Np*Nfocal, 1);
    else
        detspect{ii} = repmat(detspect{ii}, Nfocal, 1);
    end
end

% memory limit
% maxview = floor(Mlimit*2^27/(Np*Nsample*Nfocal)) * Nfocal;
maxview = 1;
Nlimit = ceil(Nview/maxview);

% ini P & Eeff
P = cell(1, Nw);
Eeff = cell(1, Nw);
switch lower(method)
    case {'default', 1, 'photoncount', 2}
        P(:) = {zeros(Np, Nview)};
        Eeff(:) = {zeros(Np, Nview)};
    case {'energyvector', 3}
        P(:) = {zeros(Np*Nview, Nsample)};
        % No Eeff
    otherwise
        % error
        error(['Unknown projection method: ' method]);
end
Pair = cell(1, Nw);

% projection on bowtie and filter in collimation
[Dmu_air, L] = flewoverbowtie(focalposition, detector.position, bowtie, filter, samplekeV);
% distance curse
distscale = detector.pixelarea(:)./(L.^2.*(pi*4));
% energy based Posibility of air
for ii = 1:Nw
    switch lower(method)
        case {'default', 1}
            % ernergy integration
            Pair{ii} = (exp(-Dmu_air).*detspect{ii}) * samplekeV';
            Pair{ii} = Pair{ii}.*distscale(:);
        case {'photoncount', 2}
            % photon counting
            Pair{ii} = sum(exp(-Dmu_air).*detspect{ii}, 2);
            Pair{ii} = Pair{ii}.*distscale(:);
        case {'energyvector', 3}
            % maintain the components on energy
            Pair{ii} = exp(-Dmu_air).*detspect{ii};
            Pair{ii} = Pair{ii}.*distscale(:);
        otherwise
            % error
            error(['Unknown projection method: ' method]);
    end
    Pair{ii}(isnan(Pair{ii})) = 0;
end
toc;

% % projection on objects
tic;
[D, mu] = projectinphantom(focalposition, detector.position, phantom, samplekeV, viewangle, couch, gantrytilt);
toc;

tic;

D = gpuArray(single(reshape(D, Np, Nview, [])));
mu = gpuArray(single(mu));
samplekeV_GPU = gpuArray(single(samplekeV));
detspect = cellfun(@(x) gpuArray(single(x)), detspect, 'UniformOutput', false);
Nlimit = gpuArray(single(Nlimit));
Np = gpuArray(single(Np));
Nfocal = gpuArray(single(Nfocal));
distscale = gpuArray(single(distscale));
maxview = gpuArray(single(maxview));

Dmu0 = gpuArray(Dmu_air);
% for i_lim = 1:Nlimit
for iview = 1:Nview
    % echo '.'
%     if echo_onoff, fprintf('.'); end
%     % view number for each step
%     if i_lim < Nlimit
%         Nview_lim = maxview;
%     else
%         Nview_lim = Nview - maxview*(Nlimit-1);
%     end
%     % index of view angles 
%     index_lim = (1:Nview_lim) + maxview*(i_lim-1);
%     index_D = (1:Np)' + (index_lim-1).*Np;
    
    % projection on objects    
    if ~isempty(D)
        Pmu = Dmu0 + squeeze(D(:, iview, :))*mu;
    else
        Pmu = Dmu0;
    end
       
    % energy based Posibility
    for ii = 1:Nw
        switch lower(method)
            case {'default', 1}
                % ernergy integration
                Pmu = exp(-Pmu).*detspect{ii};
                % for quanmtum noise
                Eeff{ii}(:, iview) = gather(reshape(sqrt((Pmu * (samplekeV'.^2))./sum(Pmu, 2)), Np, 1));
                % Pmu
                Pmu =  Pmu * samplekeV';
                Pmu = reshape(Pmu, Np*Nfocal, Nview_lim/Nfocal).*distscale(:);
                P{ii}(:, iview) = gather(reshape(Pmu, Np, 1));
%                 P{ii}(:, index_lim) = gather(reshape(reshape(Pmu*samplekeV', Np*Nfocal, Nview_lim/Nfocal).*distscale(:), Np, Nview_lim));
            case {'photoncount', 2}
                % photon counting
                Pmu = exp(-Dmu).*repmat(detspect{ii}, Nview_lim, 1);
                Pmu = sum(Pmu, 2);
                Pmu = reshape(Pmu, Np*Nfocal, Nview/Nfocal).*distscale(:);
                P{ii}(:, index_lim) = reshape(Pmu, Np, Nview_lim);
            case {'energyvector', 3}
                % maintain the components on energy
                Pmu = exp(-Dmu).*repmat(detspect{ii}, Nview_lim, 1);
                Pmu = reshape(Pmu, Np*Nfocal, []).*distscale(:);
                index_p = (1:Nview_lim*Np) + maxview*Np*(i_lim-1);
                P{ii}(index_p, :) = reshape(Pmu, Np*Nview_lim, []);
            otherwise
                % error
                error(['Unknown projection method: ' method]);
        end
    end
end
toc;

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

% end
