% run after bnc3in1cali_script1

%% step2. non-linear #2

% inputs
% BHcalitable, output of step1
% configure file
calioutputpath = 'E:\data\calibration\bh\';
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_cali.xml';
protocol_cfgfile = 'E:\matlab\CTsimulation\cali\calixml\protocol_beamharden.xml';
% rawdata
rawdata_file = struct();
rawdata_file.body = {[], [], [], []};
rawdata_file.head = {[], [], [], []};
% my data
myrawpath = 'E:\data\rawdata\bhtest\';
air_body = [myrawpath 'rawdata_air_120KV300mA_large_v1.0.raw'];
water_body = {'rawdata_water200c_120KV300mA_large_v1.0.raw', 'rawdata_water200off100_120KV300mA_large_v1.0.raw', ...
              'rawdata_water300c_120KV300mA_large_v1.0.raw', 'rawdata_water300off100_120KV300mA_large_v1.0.raw'};
rawdata_file.body{3} = {air_body, [myrawpath water_body{1}], [myrawpath water_body{2}], [myrawpath water_body{3}], ...
                        [myrawpath water_body{4}]};
% phantom
% I know the phantoms to scan are air, small water center/off and big water center/off
phantoms = {'phantom_air.xml', 'phantom_shellwater200_center.xml', 'phantom_shellwater200_off90.xml', ...
            'phantom_shellwater300_center.xml', 'phantom_shellwater300_off90.xml'};
phatompath = 'E:\matlab\CTsimulation\system\mod\phantom\';
Nphantom = length(phantoms);    % =5

% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);
Nseries = configure.protocol.seriesnumber;
% SYS
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% KV 
KV = configure.protocol.series{1}.KV;
Nw = length(KV);

% bad channel
badchannelindex = [2919 12609];
% pipe for air calibration
pipe_air = struct();
pipe_air.Log2 = struct();
pipe_air.Badchannel = struct();
pipe_air.Badchannel.badindex = badchannelindex;
pipe_air.Aircali = struct();
% pipe for noneliear correction
pipe_nl = struct();
pipe_nl.Log2 = struct();
pipe_nl.Air = struct();
pipe_nl.Badchannel.badindex = badchannelindex;
pipe_nl.Beamharden = struct();
pipe_nl.Housefield.HCscale = 1000;
pipe_nl.Databackup.dataflow = 'rawdata';
pipe_nl.Databackup.prmflow = 'recon';
pipe_nl.Axialrebin.QDO = 0;
pipe_nl.Watergoback.filter.name = 'hann';
pipe_nl.Watergoback.filter.freqscale = 1.2;
pipe_nl.Inverserebin = struct();

% ini
reconxml = struct();
reconxml.body = cell(size(phantoms));
reconxml.head = cell(size(phantoms));
% to scan data from real CT or from simulation
for iph = 1:Nphantom
    % set phantom
    phantomcfg = readcfgfile(fullfile(phatompath, phantoms{iph}));
    SYS.phantom = phantomconfigure(phantomcfg);
    % loop the series (body and head bowtie)
    for i_series = 1:Nseries
        SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
        SYS.protocol.series_index = i_series;
        % load protocol (to SYS)
        SYS = loadprotocol(SYS);
        % the bowtie is
        bowtie = lower(SYS.protocol.bowtie);
        % reconxml
        reconxml.(bowtie){iph} = reconxmloutput(SYS, 0);
        % If you want to scan on real CT, use the reconxml{iw}.protocol to scan the phantom 'phantoms{iph}'.
        % If you want to do simulation, run a simulation script like CTsimulation.m
        % But I know we have prepared the data:
        for iw = 1:Nw
            if ~isempty(rawdata_file.(bowtie){iw})
                reconxml.(bowtie){iph}{iw}.rawdata = rawdata_file.(bowtie){iw}{iph};
            else
                reconxml.(bowtie){iph}{iw}.rawdata = '';
            end
        end
        % replace pipe
        if iph==1
            % I know the 1st phantom is air
            reconxml.(bowtie){iph}.pipe = pipe_air;
            
        else
            % I know iph>1 are the water phantoms
            reconxml.(bowtie){iph}.pipe = pipe_nl;
        end
        
    end
end