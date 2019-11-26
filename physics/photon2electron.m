function Dataflow = photon2electron(SYS, Dataflow)
% photon energy to intensity

% to use
% tube
Nw = SYS.source.Wnumber;
KV = SYS.source.KV;
mA = SYS.source.mA;
mA_air = SYS.source.mA_air;
% DCB
T = SYS.datacollector.integrationtime;
gain = SYS.datacollector.DBBgain/SYS.world.refrencekeV/SYS.detector.mergescale;
Z0 = SYS.datacollector.DBBzero;
Tscale = 1000/SYS.datacollector.inttimeclock;

% constant
electric_charge = 1.602e-19;

% loop Nw
for iw = 1:Nw
    W = KV{iw}*mA{iw};
    PEscale = (T*1e-6*W/electric_charge/1000).*gain;
    % Intensity
    Dataflow.P{iw} = detectorslicemerge(Dataflow.P{iw}, SYS.detector, 'sum').*PEscale + Z0;
    % air main
    Dataflow.Pair{iw} = log2(detectorslicemerge(Dataflow.Pair{iw}, SYS.detector, 'sum') ...
        .*PEscale.*(mA_air{iw}/mA{iw})) - log2(T*Tscale);
end

end
