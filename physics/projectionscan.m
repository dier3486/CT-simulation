function Dataflow = projectionscan(SYS, method)
% the projection simulation

if nargin<2
    method = 'default';
end

% system components
source = SYS.source;
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
phantom = SYS.phantom;

% parameters of the system
focalposition = source.focalposition;
Nfocal = source.focalnumber;
Npixel = detector.Npixel;
Nslice = detector.Nslice;
Np = Npixel * Nslice;
Nw = source.Wnumber;
Mlimit = SYS.simulation.memorylimit;

% prepare the samplekeV, viewangle and couch
[samplekeV, viewangle, couch, shotindex] = scanprepare(SYS);
Nsample = length(samplekeV(:));
Nview = length(viewangle(:));

% spectrums normalize
sourcespect = SYS.source.spectrum;
for ii = 1:Nw
    sourcespect{ii} = sourcespect{ii}./sum(sourcespect{ii}.*samplekeV);
end
% detector response
detspect = cell(1, Nw);
for ii = 1:Nw
    detspect{ii} = sourcespect{ii}.*detector.spectresponse;
end
% NOTE: noly one reponse curve supported yet

% memory limit
maxview = floor(Mlimit*2^27/(Np*Nsample*Nfocal)) * Nfocal;
Nlimit = ceil(Nview/maxview);

% ini P
P = cell(1, Nw);
switch lower(method)
    case {'default', 1, 'photoncount', 2}
        P(:) = {zeros(Np, Nview)};
    case {'energyvector', 3}
        P(:) = {zeros(Np*Nview, Nsample)};
    otherwise
        % error
        error(['Unknown projection method: ' method]);
end
Pair = cell(1, Nw);

% projection on bowtie and filter in collimation
[Dmu_air, L] = flewoverbowtie(focalposition, detector.position, bowtie, filter, samplekeV);
% distance curse
distscale = detector.pixelarea./(L.^2.*(pi*4));
% energy based Posibility of air
for ii = 1:Nw
    switch lower(method)
        case {'default', 1}
            % ernergy integration
            Pair{ii} = (exp(-Dmu_air).*detspect{ii}) * samplekeV';
            Pair{ii} = Pair{ii}.*distscale(:);
        case {'photoncount', 2}
            % photon counting
            Pair{ii} = exp(-Dmu_air) * detspect{ii}';
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
    
for i_lim = 1:Nlimit
    fprintf('.')
    % view number for each step
    if i_lim < Nlimit
        Nview_lim = maxview;
    else
        Nview_lim = Nview - maxview*(Nlimit-1);
    end
    % index of view angles 
    index_lim = (1:Nview_lim) + maxview*(i_lim-1);
    % ini Dmu
    Dmu = repmat(Dmu_air, Nview_lim/Nfocal, 1);
    
    % projection on objects
    Dmu = Dmu + projectinphantom(focalposition, detector.position, phantom, samplekeV, viewangle(index_lim), couch(index_lim, :));
    
    % effective quantum patical number normalize
    % TBC
    
    % energy based Posibility
    for ii = 1:Nw
        switch lower(method)
            case {'default', 1}
                % ernergy integration
                Dmu = (exp(-Dmu).*detspect{ii}) * samplekeV';
                Dmu = reshape(Dmu, Np*Nfocal, Nview_lim/Nfocal).*distscale(:);
                P{ii}(:, index_lim) = reshape(Dmu, Np, Nview_lim);
            case {'photoncount', 2}
                % photon counting
                Dmu = exp(-Dmu) * detspect{ii}';
                Dmu = reshape(Dmu, Np*Nfocal, Nview/Nfocal).*distscale(:);
                P{ii}(:, index_lim) = reshape(Dmu, Np, Nview_lim);
            case {'energyvector', 3}
                % maintain the components on energy
                Dmu = exp(-Dmu).*detspect{ii};
                Dmu = reshape(Dmu, Np*Nfocal, []).*distscale(:);
                index_p = (1:Nview_lim*Np) + maxview*Np*(i_lim-1);
                P{ii}(index_p, :) = reshape(Dmu, Np*Nview_lim, []);
            otherwise
                % error
                error(['Unknown projection method: ' method]);
        end
        % remove nan (could due to negative thick filter)
        P{ii}(isnan(P{ii})) = 0;
    end
end

% return
Dataflow = struct();
Dataflow.samplekeV = samplekeV;
Dataflow.viewangle = viewangle;
Dataflow.couch = couch;
Dataflow.shotindex = shotindex;
Dataflow.P = P;
Dataflow.Pair = Pair;

end
