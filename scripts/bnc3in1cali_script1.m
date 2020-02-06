%% step 0. blade experiment
% to get the response curve, e.g. E:\matlab\CT\SINO\TM\detector\response_1219.mat
% skip

%% step1. beam harden #1

% inputs
% system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_basic.xml';
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_arti.xml';
protocol_cfgfile = 'E:\matlab\CTsimulation\cali\calixml\protocol_beamharden.xml';
rawdata_file = struct();
rawdata_file.empty = {[], [], 'E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_empty_v1.0.raw', []};
rawdata_file.body = {[], [], 'E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_large_v1.0.raw', []};
rawdata_file.head = {[], [], [], []};
% the rawdata_file should includes:
% {air_80kV_empty, air_100kV_empty, air_120kV_empty, air_140kV_empty}, 
% {air_80kV_body, air_100kV_body, air_120kV_body, air_140kV_body}, 
% {air_80kV_head, air_100kV_head, air_120kV_head, air_140kV_head}
% for each BH table (8 tables).
response_file = 'E:\matlab\CT\SINO\TM\detector\response_1219.mat';

% view skip in meaning the rawdata
viewskip = 200;
% bad channel
badchannelindex = [2919 12609];
% pipe
pipe_bh = struct();
pipe_bh.Log2 = [];
pipe_bh.Badchannel = struct();
pipe_bh.Badchannel.badindex = badchannelindex;

% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);
% replace the response
configure.system.detector.spectresponse = response_file;
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% The projection of air (for bowtie), simulation and real scan
Nseries = configure.protocol.seriesnumber;
P = struct();
rawmean= struct();
% loop the series
for i_series = 1:Nseries
    % I know the series 1,2,3 should be empty, body and head bowtie
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;    
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    % the bowtie is
    bowtie = lower(SYS.protocol.bowtie);
    % projection
    Data = projectionscan(SYS, 'energyvector', 0);
    % slice merge
    Nw = SYS.source.Wnumber;
    for iw = 1:Nw
        Data.Pair{iw} = ...
            detectorslicemerge(Data.Pair{iw}, SYS.detector.Npixel, SYS.detector.Nslice, SYS.detector.slicemerge, 'sum');
    end
    
    % get simulation data (ideal air)
    P.(bowtie) = Data.Pair;
    
    % get experiment data (scan air)
    reconxml = reconxmloutput(SYS, 0);
    % If you want to scan on real CT, use the reconxml{iw}.protocol to scan.
    % If you want to scan by simulation, take a look on CTsimulation.m, you should modify the SYS to employ the artifacts 
    % like quantumn noise, non-linear effects, cross-talk and so on, then run a script as in CTsimulation.m.
    % And, if you have prepared the data, do this
    rawmean.(bowtie) = cell(1, Nw);
    for iw = 1:Nw
        if ~isfield(rawdata_file, bowtie) || isempty(rawdata_file.(bowtie){iw})
            continue;
        end
        % reconcfg
        status.reconcfg = reconxml{iw};
        status.series_index = 1;
        % replace rawdata file
        status.reconcfg.rawdata = rawdata_file.(bowtie){iw};
        % replace pipe
        status.reconcfg.pipe = pipe_bh;
        % 'recon' access (loadrawdata, log2, badchannel)
        dataflow = recon_access(status);
        % mean
        rawmean.(bowtie){iw} = mean(dataflow.rawdata(:, viewskip+1:end), 2);
    end
end

% fix the bowtie thickness by fitting the simulation with experiment data
% to get the beam harden correction
BHcalitable = struct();
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
    
    % paramters to use
    Npixel = SYS.detector.Npixel;
    Nslice = max(SYS.detector.slicemerge);      % slice number after merge
    Nps = Npixel*Nslice;
    KV = SYS.protocol.KV;
    Nw = SYS.source.Wnumber;
    % material of the bowtie to fix
    mu_1 = SYS.collimation.bowtie{1}.material.mu_total(:);
    samplekeV = SYS.world.samplekeV(:);
    
    % ini the results
    BHcalitable.(bowtie) = cell(1, Nw);
    % loop KV
    for iw = 1:Nw
        if isempty(rawmean.(bowtie){iw})
            continue;
        end
        % the simulated effective empty bowtie
        Dempty = log(P.empty{iw}*samplekeV);
        % the experiment effective bowtie thickness
        Dexp = (rawmean.(bowtie){iw} - rawmean.empty{iw}).*log(2);
        % try to fit the thickness fix 'dfit' to satisfy Dbowtie(dfit)-Dempty = Dexp.
        dfit = zeros(Npixel, Nslice);
        for ipixel = 1:Nps
            Pbow_ip = P.(bowtie){iw}(ipixel, :);
            Dexp_ip = Dexp(ipixel)-Dempty(ipixel);
            if isfinite(Dexp_ip)
                dfit(ipixel) = fzero(@(x) -log(Pbow_ip*(exp(-x.*mu_1).*samplekeV))-Dexp_ip, 0);
            end
        end
        % smooth
        for islice = 1:Nslice
            dfit(:, islice) = smooth(dfit(:,islice), 0.05, 'rloess');
        end
        % set KV
        SYS.protocol.KV = KV(iw);
        % reload protocol
        SYS = loadprotocol(SYS);
        % add effect filter
        Nfilt = length(SYS.collimation.filter);
        SYS.collimation.filter{Nfilt+1} = struct();
        SYS.collimation.filter{Nfilt+1}.effect = true;
        SYS.collimation.filter{Nfilt+1}.thickness = dfit(:);
        SYS.collimation.filter{Nfilt+1}.material = SYS.collimation.bowtie{1}.material;
        
        % BH cali
        BHcorr = simuBHcali(SYS, 4);
        BHcalitable.(bowtie){iw} = BHcorr{1};
        % save table
        corrfile = fullfile(SYS.output.path, [SYS.output.files.beamharden{1} '.corr']);
        cfgfile = cfgmatchrule(corrfile, SYS.path.IOstandard);
        corrcfg = readcfgfile(cfgfile);
        packstruct(BHcorr{1}, corrcfg, corrfile);
    end
end
