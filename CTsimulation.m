function reconxml = CTsimulation(configure_file)
% main function of the CT simulation
% CTsimulation(configure_file);
% INPUT:
%   configure_file      the configure .xml file 
% You may find a sample of the configure file in ~\system\mod\sample_configure.xml
% OUTPUT: output files in the output path set in configure file

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% where am I
mainfile = which('CTsimulation');
rootpath = fileparts(mainfile);
% add path
addpath(genpath(rootpath));

% load configure file
if ischar(configure_file) || isstring(configure_file)
    configure = readcfgfile(configure_file);
elseif isstruct(configure_file)
    configure = configure_file;
else
    error('Illegal input in CTsimulation.');
end
configure = configureclean(configure);

% system configure
fprintf('system configure...');
SYS = systemconfigure(configure.system);
% phantom configure
if isfield(configure, 'phantom')
    SYS.phantom = phantomconfigure(configure.phantom);
end
% simulation prepare (load materials)
SYS = systemprepare(SYS);
fprintf(' done\n');

% serie number
Nseries = configure.protocol.seriesnumber;
% ini return
reconxml = cell(1, Nseries);
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
    % no scatter now
    1;
    % to intensity
    fprintf('  to intensity...');
    Data = photon2electron(SYS, Data);
    fprintf(' done\n');
    % output rawdata, corr table and recon xml
    fprintf('  output to datapath...');
    reconxml{i_series} = simuresultsoutput(SYS, Data);
    fprintf(' done.\n');
end

end