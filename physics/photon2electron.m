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
gain = SYS.datacollector.DBBgain/SYS.detector.mergescale;
Z0 = SYS.datacollector.DBBzero;
Tscale = 1000/SYS.datacollector.inttimeclock;

% constant
electric_charge = 1.602e-19;

% loop Nw
for iw = 1:Nw
    W = KV{iw}*mA{iw};
    PEscale = (T*1e-6*W/electric_charge/1000).*gain;
    % Intensity
    Dataflow.P{iw} = DCBslicemerge(Dataflow.P{iw}, SYS.detector).*PEscale + Z0;
    % air main
    Dataflow.Pair{iw} = log2(DCBslicemerge(Dataflow.Pair{iw}, SYS.detector) ...
        .*PEscale.*(mA_air{iw}/mA{iw})) - log2(T*Tscale);
end

end


function Pout = DCBslicemerge(Pin, detector)
% merge the slices

Nslice = detector.Nslice;
Npixel = detector.Npixel;

if all(detector.slicemerge == 1:Nslice)
    % skip
    Pout = Pin;
    return;
end

% to merge
Nmergedslice = max(detector.slicemerge);
Pin = reshape(Pin, Npixel, Nslice, []);
Pout = zeros(Npixel, Nmergedslice, size(Pin, 3));
for ii = 1:Nslice
    index_ii = detector.slicemerge(ii);
    Pout(:, index_ii, :) = Pout(:, index_ii, :) + Pin(:, ii, :);
end
Pout = reshape(Pout, Npixel*Nmergedslice, []);

end