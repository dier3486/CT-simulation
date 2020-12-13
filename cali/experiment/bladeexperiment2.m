% Cu blade experiment for fitting detector response

% addpth
CTsimupath = '../../';
addpath(genpath(CTsimupath));

% CT system
% configure_file = 'E:\matlab\calibration\system\configure_cali.xml';
configure_file = 'E:\matlab\CT\SINO\PG\configure_PG_cali.xml';
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

% I know the 1st series is empty bowtie
i_series = 1;

SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
SYS.protocol.seriesindex = i_series;
% load protocol (to SYS)
SYS = loadprotocol(SYS);

% projection 
Data = projectionscan(SYS, 'energyvector');

% I know the blades are
bladestep = 0.5;
blademaxthick = 20; 
blades = bladestep: bladestep: blademaxthick;

Nbld = length(blades);
Nw = SYS.source.Wnumber;
Npixel = double(SYS.detector.Npixel);
Nslice = double(SYS.detector.Nslice);
Np = Npixel*Nslice;
Nsmp = length(SYS.world.samplekeV);

% get the simulation results 
Pbld = cell(1, Nw);
for iw = 1:Nw
    Pbld{iw} = ((Data.P{iw}(:)./Data.Pair{iw}(:)).^blades).*Data.Pair{iw}(:);
    Pbld{iw}(isnan(Pbld{iw})) = 0;
    Pbld{iw} = reshape(Pbld{iw}, Np, Nsmp, Nbld);
end

% samplekeV & detangle
samplekeV = SYS.world.samplekeV;
[detangle, eqangle] = detpos2angle(SYS.detector.position(1:Npixel, :), SYS.source.focalposition(1,:));

% load experiment data
% experimentdata = 'F:\data-Dier.Z\Cublades\Data1218.mat';
experimentdata = 'F:\data-Dier.Z\PG\Cu\Data0216.mat';
blddata = load(experimentdata);

expbld = cell(1, Nw);
expbld(:) = {zeros(Np, Nbld)};
expair = cell(1, Nw);
expair(:) = {zeros(Np, 1)};
Ndata = length(blddata.Data);
for ii = 1:Ndata
    % KV
    switch blddata.Data(ii).KV
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
    % obj
    obj = regexp(blddata.Data(ii).object, '[a-z_A-Z]+', 'match');
    obj = obj{1};
    switch lower(obj)
        case 'air'
            % air data
            expair{windex} = blddata.Data(ii).rawmean(:) - blddata.Data(ii).offsetmean(:);
        case 'cu'
            if strcmpi(blddata.Data(ii).focal, 'small')
                thick = regexp(blddata.Data(ii).object, '[^a-z_A-Z]+', 'match');
                thick = str2double(thick{1});
                thick_index = round(thick/bladestep);
                expbld{windex}(:, thick_index) = blddata.Data(ii).rawmean(:) - blddata.Data(ii).offsetmean(:);
            end
    end
    %
end

% expdata_bow = 'F:\data-Dier.Z\stepdata\Data1205.mat';
% bowdata = load(expdata_bow);
% Ndata = length(bowdata.Data);
% expbow = cell(1, Nw);
% expbow(:) = {zeros(Np, 3)};
% for ii = 1:Ndata
%     % KV
%     switch bowdata.Data(ii).KV
%         case 80
%             windex = 1;
%         case 100
%             windex = 2;
%         case 120
%             windex = 3;
%         case 140
%             windex = 4;
%         otherwise
%             windex = [];
%     end
%     % bowtie
%     switch lower(bowdata.Data(ii).bowtie)
%         case 'empty'
%             bindex = 0;
%         case 'body'
%             bindex = 1;
%         case 'head'
%             bindex = 2;
%     end
%     % obj
%     obj = regexp(bowdata.Data(ii).object, '[a-z_A-Z]+', 'match');
%     obj = obj{1};
%     switch obj
%         case 'air'
%             % air data
%             expbow{windex}(:, bindex+1) = bowdata.Data(ii).rawmean(:) - bowdata.Data(ii).offsetmean(:);
%         otherwise
%             1;
%     end
%     
% end





