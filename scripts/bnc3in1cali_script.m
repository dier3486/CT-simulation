%% step 0. blade experiment
% to get the response curve, e.g. E:\matlab\CT\SINO\TM\detector\response_1219.mat
% skip

%% step1. beam harden #1

% inputs
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_basic.xml';
protocol_cfgfile = 'E:\matlab\calibration\system\protocol_beamharden.xml';
rawdata_file.empty = {[], [], 'E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_empty_v1.0.raw', []};
rawdata_file.body = {[], [], 'E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_large_v1.0.raw', []};
rawdata_file.head = {[], [], [], []};
% the rawdata_file should includes:
% {air_80kV_empty, air_100kV_empty, air_120kV_empty, air_140kV_empty}, 
% {air_80kV_body, air_100kV_body, air_120kV_body, air_140kV_body}, 
% {air_80kV_head, air_100kV_head, air_120kV_head, air_140kV_head}
% for each BH table (8 tables).
response_file = 'E:\matlab\CT\SINO\TM\detector\response_1219.mat';

% % load experiment data (or do scan to get these data, or do simulation to get these data)
% rawmean= struct();
% viewskip = 100;
% smallcfg.IOstandard = [];
% rawfields = fieldnames(rawdata_file);
% for ifield = 1:length(rawfields)
%     field_ii = rawfields{ifield};
%     Nw = length(rawdata_file.(field_ii));
%     rawmean.(field_ii) = cell(1, Nw);
%     for iw = 1:Nw
%         if ~isempty(rawdata_file.(field_ii){iw})
%             smallcfg.rawdata = rawdata_file.(field_ii){iw};
%             % quick data reading
%             dataflow = readrawdata(smallcfg);
%             % log2
%             dataflow = reconnode_log2(dataflow, [], []);
%             % mean
%             rawmean.(field_ii){iw} = mean(dataflow.rawdata(:, viewskip+1:end), 2);
%         end
%     end
% end

% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);
% replace the response
configure.system.detector.reponse = response_file;
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% simulate the projection of air (for bowtie)
Nseries = configure.protocol.seriesnumber;
P = struct();
% loop the series
for i_series = 1:Nseries
    % I know the series 1,2,3 should be empty, body and head bowtie
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;    
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);  
    % projection
    Data = projectionscan(SYS, 'energyvector', 0);
    % slice merge
    Nw = SYS.source.Wnumber;
    for iw = 1:Nw
        Data.Pair{iw} = ...
            detectorslicemerge(Data.Pair{iw}, SYS.detector.Npixel, SYS.detector.Nslice, SYS.detector.slicemerge, 'sum');
    end
    
    % get simulation data (ideal air)
    bowtie = lower(SYS.protocol.bowtie);
    P.(bowtie) = Data.Pair;
    
    % get experiment data
    bowtie = lower(SYS.protocol.bowtie);
    reconxml = reconxmloutput(SYS, false);
    % If you want to scan on real CT, use the reconxml{iw}.protocol to scan.
    % If you want to scan by simulation, take a look on CTsimulation.m, you should modify the SYS to employ the artifacts 
    % like quantumn noise, non-linear effects, cross-talk and so on, then run a script as in CTsimulation.m.
    % And, if you have prepared the data, do this
    for iw = 1:Nw
        if ~isfield(rawdata_file, bowtie) || isempty(rawdata_file.(bowtie){iw})
            continue;
        end
        % reconcfg
        status.reconcfg = reconxml{iw};
        status.series_index = 1;
        % fix pipe
        
        % initial
        dataflow = struct();
        prmflow = struct();
        
        
    end
end

% samplekeV & detangle
% samplekeV = SYS.world.samplekeV;
% Npixel = SYS.detector.Npixel;
% [detangle, eqangle] = detpos2angle(SYS.detector.position(1:Npixel, :), SYS.source.focalposition(1,:));

% fix the bowtie thickness by fitting the simulation with experiment data
for i_series = 1:Nseries
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    bowtie = lower(SYS.protocol.bowtie);
    if strcmpi(bowtie, 'empty')
        % skip the BH correction for empty bowtie
        continue;
    end
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    % material of the bowtie to fix
    mu_1 = SYS.collimation.bowtie{1}.material.mu_total(:);
    samplekeV = SYS.world.samplekeV(:);
    % paramters to use
    Npixel = SYS.detector.Npixel;
    Nslice = max(SYS.detector.slicemerge);      % slice number after merge
    Nps = Npixel*Nslice;
    
    % ini the fit results
    dfit = cell(1, Nw);
%     dfit(:) = {zeros(Nps, 1)};
    % loop KV
    Nw = SYS.source.Wnumber;
    for iw = 1:Nw
        if isempty(rawmean.(bowtie){iw})
            continue;
        end
        % the simulated effective empty bowtie
        Dempty = log(P.empty{iw}*samplekeV);
        % the experiment effective bowtie thickness
        Dexp = (rawmean.(bowtie){iw} - rawmean.empty{iw}).*log(2);
        % try to fit the thickness-fix x to satisfy Dbowtie(x)-Dempty = Dexp.
        dfit{iw} = zeros(Nps, 1);
        for ipixel = 1:Nps
            Pbow_ip = P.(bowtie){iw}(ipixel, :);
            Dexp_ip = Dexp(ipixel)-Dempty(ipixel);
            if isfinite(Dexp_ip)
                dfit{iw}(ipixel) = fzero(@(x) -log(Pbow_ip*(exp(-x.*mu_1).*samplekeV))-Dexp_ip, 0);
            end
        end
        1;
    end
    
end


%% step2.
