clear;
addpath(genpath('../'));

% configure_file = '..\system\mod\sample_configure.xml';
configure_file = 'E:\matlab\CT\SINO\TM\configure_1.xml';

% load configure file
configure = readcfgfile(configure_file);
configure = configureclean(configure);

% system configure
SYS = systemconfigure(configure.system);
% phantom configure
SYS.phantom = phantomconfigure(configure.phantom);
% simulation prepare (load materials)
SYS = systemprepare(SYS);

Nseries = configure.protocol.seriesnumber;
% loop the series
for i_series = 1:Nseries
    % to play i-th series
    % protocol configure
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    % projection
    Data = projectionscan(SYS);
    % to intansity
    Data = photon2electron(SYS, Data);
    % output to rawdata
    simuresultsoutput(SYS, Data);
end
