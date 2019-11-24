function CTsimulation(configure_file)

% where am I
mainfile = which('CTsimulation');
rootpath = fileparts(mainfile);
% add path
addpath(genpath(rootpath));

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
    % load protocol
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    SYS.protocol.series_index = i_series;
    SYS = loadprotocol(SYS);
    % projection
    [P, Pair] = projectionscan(SYS);
    %
end



end