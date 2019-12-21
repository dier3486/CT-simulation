% Cu blade experiment for fitting detector response

% add path
CTsimupath = '../../';
addpath(genpath(CTsimupath));

% CT system
configure_file = 'E:\matlab\calibration\system\configure_cali.xml';
% read configure file
configure = readcfgfile(configure_file);
% load configure
configure = configureclean(configure);
% system configure
SYS = systemconfigure(configure.system);
% phantom configure
SYS.phantom = phantomconfigure(configure.phantom);
% simulation prepare (load materials)
SYS = systemprepare(SYS);

% simulation of the bowties
% loop the series (for each bowtie)
Nseries = configure.protocol.seriesnumber;
Nbow = Nseries;
for i_series = 1:Nseries
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    % ini 
    if i_series == 1
        Nw = SYS.source.Wnumber;
        Npixel = double(SYS.detector.Npixel);
        Nslice = double(SYS.detector.Nslice);
        Np = Npixel*Nslice;
        Nsmp = length(SYS.world.samplekeV);
        Pbow = cell(1, Nw);
        Pbow(:) = {zeros(Np, Nsmp, Nbow)};
    end
    % projection
    Data = projectionscan(SYS, 'energyvector');
    
    % I know iw is KV, i_series is bowtie
    for iw = 1:Nw
        Pbow{iw}(:, :, i_series) = Data.Pair{iw};
    end
    % copy the collimation
    collim(i_series) = SYS.collimation;
end

% samplekeV & detangle
samplekeV = SYS.world.samplekeV;
[detangle, eqangle] = detpos2angle(SYS.detector.position(1:Npixel, :), SYS.source.focalposition(1,:));


% experiment
expdata_bow = 'F:\data-Dier.Z\stepdata\Data1205.mat';
bowdata = load(expdata_bow);
Ndata = length(bowdata.Data);
expbow = cell(1, Nw);
expbow(:) = {zeros(Np, 3)};
for ii = 1:Ndata
    % KV
    switch bowdata.Data(ii).KV
        case 80
            windex = 1;
        case 100
            windex = 2;
        case 120
            windex = 3;
        case 140
            windex = 4;
        otherwise
            windex = [];
    end
    % bowtie
    switch lower(bowdata.Data(ii).bowtie)
        case 'empty'
            bindex = 1;
        case {'body', 'large'}
            bindex = 2;
        case {'head', 'small'}
            bindex = 3;
    end
    % obj
    obj = regexp(bowdata.Data(ii).object, '[a-z_A-Z]+', 'match');
    obj = obj{1};
    switch obj
        case 'air'
            % air data
            expbow{windex}(:, bindex) = bowdata.Data(ii).rawmean(:) - bowdata.Data(ii).offsetmean(:);
        otherwise
            1;
    end
end

% the returns are
% Pbow, expbow, collim, detangle and eqangle