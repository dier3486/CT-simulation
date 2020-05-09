function Dataflow = airprojection(SYS, method)
% the projection simulation

if nargin<2 || isempty(method)
    method = 'default';
end

% system components
source = SYS.source;
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;

% parameters of the system
focalposition = source.focalposition;
Nfocal = source.focalnumber;
Npixel = double(detector.Npixel);
Nslice = double(detector.Nslice);
Np = Npixel * Nslice;
Nw = source.Wnumber;

% prepare the samplekeV
samplekeV = scanprepare(SYS);

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

% ini
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
    
    % slice merge (only for airprojection)
    Pair{ii} = detectorslicemerge(Pair{ii}, Npixel, Nslice, detector.slicemerge, 'sum');
end

% return
Dataflow = struct();
Dataflow.samplekeV = samplekeV;
Dataflow.Pair = Pair;

end
