function [Intensity, air_main] = photon2electron(SYS, P, Pair)
% photon energy to intensity

% to use
DCB = SYS.datacollector;
Nw = SYS.source.Wnumber;
KV = SYS.source.KV;
mA = SYS.source.mA;
mA_air = SYS.source.mA_air;
T = SYS.protocol.integrationtime;
gain = SYS.datacollector.gain;

% constant
electric_charge = 1.602e-19;

for iw = 1:Nw
    W = KV{iw}*mA{iw};
    W_air = KV{iw}*mA_air{iw};
    Pscale = (T*1e-6*W/electric_charge/1000).*gain;
end




W = KV*mA;
W_air = W;
T = SYS.protocol.integrationtime;
gain = 0.1;
Z0 = 16384;
maxanglecode = 69120;