function CTsimulation(configure_file)
% main function of the CT simulation

% where am I
mainfile = which('CTsimulation');
rootpath = fileparts(mainfile);
% add path
addpath(genpath(rootpath));

% load configure file
configure = readcfgfile(configure_file);
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
    SYS.protocol.series_index = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    fprintf(' done\n');
    % projection
    fprintf('  projection');
    Data = projectionscan(SYS);
    fprintf(' done\n');
    % no scatter now
    1;
    % no quatum noise now
    1;
    % to intensity
    fprintf('  output to datapath...');
    Data = photon2electron(SYS, Data);
    % output rawdata, corr table and recon xml
    simuresultsoutput(SYS, Data);
    fprintf(' done.\n');
end

end