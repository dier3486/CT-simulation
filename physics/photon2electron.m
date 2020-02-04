function Dataflow = photon2electron(SYS, Dataflow)
% photon energy to intensity

% to use
% subs
detector = SYS.detector;
% tube
Nw = SYS.source.Wnumber;
KV = SYS.source.KV;
mA = SYS.source.mA;
mA_air = SYS.source.mA_air;
% DCB
T = SYS.datacollector.integrationtime;
gain = SYS.datacollector.DBBgain/SYS.world.referencekeV/detector.mergescale;
Z0 = SYS.datacollector.DBBzero;
Tscale = 1000/SYS.datacollector.inttimeclock;

% constant
electric_charge = 1.602e-19;
epeffect = 0.01;

% intensity scale
for iw = 1:Nw
    W = KV{iw}*mA{iw};
    PEscale = T*1e-6*W/electric_charge/1000*epeffect;
    % I know T*1e-6 is in sec, W is in KV*mA=V*Q/sec, electric_charge*1000(V) is in KeV, epeffect is 1%,
    % therefore PEscale is in V*Q/KeV, which will scale the KeV normed P in counting number on KeV,
    % P*PEscale/Ephonton is the phonton number where Ephonton is in KeV.
    
    % Intensity
    Dataflow.P{iw} = Dataflow.P{iw}.*PEscale;
    Dataflow.Pair{iw} = Dataflow.Pair{iw}.*PEscale.*(mA_air{iw}/mA{iw});
end

% offfocal
if SYS.simulation.offfocal && isfield(SYS.source, 'offfocalintensity')
    fprintf(' Off-focal...');
    offintensity = SYS.source.offfocalintensity;
    for iw = 1:Nw
        Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity) + offfocalpseudoscan(SYS, Dataflow.P{iw});
        Dataflow.Pair{iw} = Dataflow.Pair{iw}.*(1-offintensity) + offfocalpseudoscan(SYS, Dataflow.Pair{iw});
        % or
        % Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity) + offfocalpseudoscan2(SYS, Dataflow.P{iw}, Dataflow.Pair{iw});
    end
end

% Quantum noise
if SYS.simulation.quantumnoise && ~isempty(Dataflow.Eeff{iw})
    fprintf(' Quantum noise...');
    for iw = 1:Nw
        % Dataflow.P{iw} = poissrnd(Dataflow.P{iw}./Dataflow.Eeff{iw}).*Dataflow.Eeff{iw};
        Dataflow.P{iw} = (poissrnd(Dataflow.P{iw}./Dataflow.Eeff{iw})+rand(size(Dataflow.P{iw}))-0.5).*Dataflow.Eeff{iw};
    end
    % don't put quantum noise on Pair
end

% cross talk
if SYS.simulation.crosstalk && isfield(detector, 'crossmatrix')
    fprintf(' Crosstalk...');
    for iw = 1:Nw
        Dataflow.P{iw} = detector.crossmatrix\Dataflow.P{iw};
        Dataflow.Pair{iw} = detector.crossmatrix\Dataflow.Pair{iw};
    end
end

% slice merge
for iw = 1:Nw
    Dataflow.P{iw} = detectorslicemerge(Dataflow.P{iw}, detector.Npixel, detector.Nslice, detector.slicemerge, 'sum');
    % DBB gain
    Dataflow.P{iw} = Dataflow.P{iw}.*gain + Z0;
    % air main
    Dataflow.Pair{iw} = -log2(detectorslicemerge(Dataflow.Pair{iw}, detector.Npixel, detector.Nslice, ...
                        detector.slicemerge, 'sum').*gain) + log2(T*Tscale);
end

end
