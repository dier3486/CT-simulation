function Dataflow = photon2electron(SYS, Dataflow, echo_onoff)
% photon energy to intensity

if nargin<3
    echo_onoff = true;
end

% to use
Mlimit = SYS.simulation.memorylimit;
% subs
detector = SYS.detector;
Npixel = double(detector.Npixel);
Nslice = double(detector.Nslice);
Np = Npixel * Nslice;
% tube
Nw = SYS.source.Wnumber;
Nfocal = SYS.source.focalnumber;
KV = SYS.source.KV;
mA = SYS.source.mA;
mA_air = SYS.source.mA_air;
% DCB
T = SYS.datacollector.integrationtime;
DBBgain = SYS.datacollector.DBBgain/SYS.world.referencekeV/detector.mergescale;
Quantumgain = SYS.datacollector.Quantumgain;
Z0 = SYS.datacollector.DBBzero;
Tscale = 1000/SYS.datacollector.inttimeclock;
% console
% iblock = Dataflow.iblock;

% constant
electric_charge = 1.602e-19;
epeffect = 0.01;

% intensity scale
for iw = 1:Nw
    W = KV{iw}*mA{iw};
    PEscale = T*1e-6*W/electric_charge/1000*epeffect*Quantumgain;
    % I know T*1e-6 is in sec, W is in KV*mA=V*Q/sec, electric_charge*1000(V) is in KeV, epeffect is 1%,
    % therefore PEscale is in V*Q/KeV, which will scale the KeV normed P in counting number on KeV,
    % P*PEscale/Ephonton is the phonton number where Ephonton is in KeV.
    
    % Intensity
    Dataflow.P{iw} = Dataflow.P{iw}.*PEscale;
    Dataflow.Pair{iw} = Dataflow.Pair{iw}.*PEscale.*(mA_air{iw}/mA{iw});
end

% % randon mA (tmp)
% Dataflow.mA = cell(Nw, 1);
% for iw = 1:Nw
%     mAbyview = (rand(1, size(Dataflow.P{iw}, 2))-0.5).*0.05 + 1;
%     Dataflow.P{iw} = Dataflow.P{iw}.*mAbyview;
%     Dataflow.mA{iw} = mAbyview.*mA{iw};
% end

% % offfocal (deleted)
% if SYS.simulation.offfocal && isfield(SYS.source, 'offfocalintensity')
%     % sorry, only for axial
%     if strcmpi(SYS.protocol.scan, 'Axial')
%         % echo 'Off-focal'
%         if echo_onoff, fprintf(' Off-focal...'); end
%         % + off-focal
%         offintensity = SYS.source.offfocalintensity;
%         for iw = 1:Nw
%             Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity) + offfocalpseudoscan(SYS, Dataflow.P{iw});
%             Dataflow.Pair{iw} = Dataflow.Pair{iw}.*(1-offintensity) + offfocalpseudoscan(SYS, Dataflow.Pair{iw});
%             % or
%             % Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity) + offfocalpseudoscan2(SYS, Dataflow.P{iw}, Dataflow.Pair{iw});
%         end
%     end
% end

% Quantum noise
if SYS.simulation.quantumnoise && ~isempty(Dataflow.Eeff{iw})
    % echo 'Quantum noise'
    if echo_onoff, fprintf(' Quantum noise...'); end
    % memory limit
    maxview = floor(Mlimit*2^24/Np/2);
    % loop kVmA
    for iw = 1:Nw
        Nview = size(Dataflow.P{iw}, 2);
        Nlimit = ceil(Nview/maxview);
        for i_lim = 1:Nlimit
            % echo '.'
            if echo_onoff, fprintf('.'); end
            % view number for each step
            if i_lim < Nlimit
                Nview_lim = maxview;
            else
                Nview_lim = Nview - maxview*(Nlimit-1);
            end
            % index of view angles 
            index_lim = (1:Nview_lim) + maxview*(i_lim-1);
            % Dataflow.P{iw} = poissrnd(Dataflow.P{iw}(:, index_lim)./Dataflow.Eeff{iw}).*Dataflow.Eeff{iw};
            Dataflow.P{iw}(:, index_lim) = (poissrnd(Dataflow.P{iw}(:, index_lim)./Dataflow.Eeff{iw}(:, index_lim)) + ...
                rand(size(Dataflow.P{iw}(:, index_lim)))-0.5).*Dataflow.Eeff{iw}(:, index_lim);
        end
    end
    % don't put quantum noise on Pair
end

% cross talk
if SYS.simulation.crosstalk && isfield(detector, 'crossmatrix')
    % echo 'Crosstalk'
    if echo_onoff, fprintf(' Crosstalk...'); end
    for iw = 1:Nw
        Dataflow.P{iw} = detector.crossmatrix\Dataflow.P{iw};
        Dataflow.Pair{iw} = detector.crossmatrix\Dataflow.Pair{iw};
    end
end

% slice merge
for iw = 1:Nw
    Dataflow.P{iw} = detectorslicemerge(Dataflow.P{iw}, detector.Npixel, detector.Nslice, detector.slicemerge, 'sum');
    % DBB gain
    Dataflow.P{iw} = Dataflow.P{iw}.*DBBgain + Z0;
    % air main
    Dataflow.Pair{iw} = -log2(detectorslicemerge(Dataflow.Pair{iw}, detector.Npixel, detector.Nslice, ...
                        detector.slicemerge, 'sum').*DBBgain) + log2(T*Tscale);
end

end
