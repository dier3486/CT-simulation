%% step 0. blade experiment
% to get the response curve, e.g. E:\matlab\CT\SINO\TM\detector\response_1219.mat
% skip

%% step1. beam harden #1

% inputs
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_basic.xml';
protocol_cfgfile = 'E:\matlab\calibration\system\protocol_beamharden.xml';
rawdata_file.empty = {[], [], 'rawdata_staticair_120KV200mA_empty_v1.0.raw', []};
rawdata_file.body = {[], [], 'rawdata_staticair_120KV200mA_large_v1.0.raw', []};
rawdata_file.head = {[], [], [], []};
% the rawdata_file should includes:
% {air_80kV_empty, air_100kV_empty, air_120kV_empty, air_140kV_empty}, 
% {air_80kV_body, air_100kV_body, air_120kV_body, air_140kV_body}, 
% {air_80kV_head, air_100kV_head, air_120kV_head, air_140kV_head}
% for each BH table (8 tables).
response_file = 'E:\matlab\CT\SINO\TM\detector\response_1219.mat';

% load experiment data
Nw = length(rawdata_file.empty);
rawdata.empty = cell(1, Nw);
rawdata.body = cell(1, Nw);
rawdata.head = cell(1, Nw);
for iw = 1:Nw
    1;
end


% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% load response
rsp = load(response_file);

% simulate the projection of air
Nseries = configure.protocol.seriesnumber;
% loop the series
for i_series = 1:Nseries
    % I know the series 1,2,3 should be empty, body and head bowtie
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    
    % projection
    Data = projectionscan(SYS, 'energyvector', 0);
    switch lower(SYS.protocol.bowtie)
        case 'empty'
            Pempty = Data.Pair;
        case 'body'
            Pbody = Data.Pair;
        case 'head'
            Phead = Data.Pair;
    end
end

% samplekeV & detangle
samplekeV = SYS.world.samplekeV;
[detangle, eqangle] = detpos2angle(SYS.detector.position(1:Npixel, :), SYS.source.focalposition(1,:));

% go back to bowtie
for i_series = 1:Nseries
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    if strcmpi(SYS.protocol.bowtie, 'empty')
        continue;
    end
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    
    % loop KV
    Nw = SYS.source.Wnumber;
    for iw = 1:Nw
        Dair = log(Pbow{iw}(:,:,1)*(r.*samplekeV(:)));
    end
    
end


%% step2.
