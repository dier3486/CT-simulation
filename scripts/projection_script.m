clear;
% addpath(genpath('../'));


% configure_file = 'D:\data\simulation\aircali_test\TY\mod\configure.xml';
configure_file = 'E:\matlab\CT\SINO\TM\configure_1.xml';
% configure_file = 'D:\matlab\ct\BCT16\BHtest\configure_BH1.xml';
% load configure file
configure = readcfgfile(configure_file);

% % configure sample
% mainfile = which('CTsimulation');
% rootpath = [fileparts(mainfile) '\'];
% configure.system = systemcfgsample(rootpath);
% configure.phantom = phantomcfgsample();
% configure.protocol = protocolcfgsample();
% % save the samples
% system_cfg.system = configure.system;
% struct2xml(system_cfg, [rootpath 'system\mod\sample_system.xml']);
% phantom_cfg.phantom = configure.phantom;
% struct2xml(phantom_cfg, [rootpath 'system\mod\sample_phantom.xml']);
% protocol_cfg.protocol = configure.protocol;
% struct2xml(protocol_cfg, [rootpath 'system\mod\sample_protocol.xml']);

% cfg clean
configure = configureclean(configure);

% system configure
fprintf('system configure...');
SYS = systemconfigure(configure.system);
% phantom configure
SYS.phantom = phantomconfigure(configure.phantom);
% simulation prepare (load materials)
SYS = systemprepare(SYS);
fprintf(' done\n');
Nseries = configure.protocol.seriesnumber;
% loop the series
for i_series = 1:Nseries
    % to play i-th series
    fprintf('play series %d\n', i_series);
    % protocol configure
    fprintf('  load protocol...');
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.seriesindex = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    fprintf(' done\n');
    % projection
    fprintf('  projection');
    Data = projectionscan(SYS);
    fprintf(' done\n');
    % to intensity
    fprintf('  to intensity...');
    Data = photon2electron(SYS, Data);
    fprintf(' done\n');
    % output to rawdata
    fprintf('  output data...');
    simuresultsoutput(SYS, Data);
    fprintf(' done.\n');
end
