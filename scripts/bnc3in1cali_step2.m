% step2. nonlinear#1, crosstalk#1

% prepare data
% the protocols to loop
toloop = struct();
toloop.bowtie = {'body', 'head'};
toloop.focalsize = {'small', 'large'};
toloop.collimator = {'32x0.625'};
toloop.KV = [80 100 120 140];

% rawdata file path
% % air
% filepath = struct();
% filepath.air.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582870883044.0_AIR';
% filepath.air.namekey = '';
% filepath.air.skiptag = 'collimator';
% water 20cm ISO
filepath.water200c.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582884571221.0_WATER22_ISO';
filepath.water200c.namekey = '';
filepath.water200c.skiptag = 'collimator';
% water 20cm off
filepath.water200off.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582876134042.0_WATER22_9CM';
filepath.water200off.namekey = '';
filepath.water200off.skiptag = 'collimator';
% water 30cm ISO
filepath.water300c.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582883779173.0_WATER32_ISO';
filepath.water300c.namekey = '';
filepath.water300c.skiptag = 'collimator';
% water 30cm off
filepath.water300off.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582877389092.0_WATER32_9CM';
filepath.water300off.namekey = '';
filepath.water300off.skiptag = 'collimator';
% file ext
fileext = '.pd';
% get file names
datafile_nl = calidataprepare(toloop, filepath, fileext);

% % calibration tables path
% corrpath.beamharden.path = 'E:\matlab\CT\SINO\PG\calibration';
% corrpath.beamharden.namekey = {'beamharden', 'bhcali'};
% corrpath.beamharden.skiptag = 'focalsize';
% % calitable ext
% corrext = '.corr';
% % get calibration tables (beamharden)
% datafile_nl = calicorrprepare(datafile_nl, corrpath, corrext);

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PG\Nonlinearcali#1_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'E:\matlab\CT\SINO\PG\calibration\';
% % namekey
% namekey = 'none#1';

% calibration paramters
% bad channel (shall be a corr table)
badchannelindex = [];
% off-focal corr (shall be a corr table)
Offfocal = struct();
Offfocal.offintensity = [0.007 0.001];
Offfocal.offwidth = [65 95];
Offfocal.offedge = [0.6 0.6];
Offfocal.ratescale = [0.8 0.8];
% water go back to get ideal water (shall be fix for each machine version)
Watergoback = struct();
Watergoback.QDO = false;
Watergoback.filter.name = 'hann';
Watergoback.filter.freqscale = 1.5;
Watergoback.span = 30;
% Watergoback.offfocal = 'deep';
% Watergoback.offfocal = 'weak';
Watergoback.offfocal = 'none';
Watergoback.offplot = true;
% nonlinear cali
nonlinearcali = struct();
nonlinearcali.corrversion = 'v1.11';
% crosstalk cali (shall be fix each machine version)
crosstalkcali = struct();
crosstalkcali.Npixelpermod = 16;
crosstalkcali.Nmerge = 4;
% crosstalkcali.istointensity = true;   (default = true)
crosstalkcali.corrversion = 'v1.11';

% debug
datafile_nl = datafile_nl(11:12);
% I know the 11 is 120KV, head bowtie, small focal

% loop the protocols
Nprotocol = size(datafile_nl(:), 1);
for ii = 1:Nprotocol
    if isempty(datafile_nl(ii).filename)
        continue;
    end
    % set the values in cali xml
    calixml = calibase;
    % set the paramters to pipe line
    for jj = 1:2
        calixml.recon{jj}.outputpath = calioutputpath;
        calixml.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml.recon{jj}.pipe.Offfocal = structmerge(Offfocal, calixml.recon{jj}.pipe.Offfocal);
        calixml.recon{jj}.pipe.Idealwater = structmerge(Watergoback, calixml.recon{jj}.pipe.Idealwater);
    end
    % 1st, water 20cm off
    calixml.recon{1}.rawdata = datafile_nl(ii).filename.water200off;
    calixml.recon{1}.protocol.namekey = 'water200off';
    % 2nd, water 30cm off
    calixml.recon{2}.rawdata = datafile_nl(ii).filename.water300off;
    calixml.recon{2}.protocol.namekey = 'water300off';
    % cali
    calixml.recon{2}.pipe.nonlinearcali = structmerge(nonlinearcali, calixml.recon{2}.pipe.nonlinearcali);
    calixml.recon{2}.pipe.crosstalkcali = structmerge(crosstalkcali, calixml.recon{2}.pipe.crosstalkcali);
    
    % echo
    fprintf('Nonlinear Calibration #1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_nl(ii).collimator, datafile_nl(ii).bowtie, datafile_nl(ii).KV, datafile_nl(ii).focalsize);
    
    % run the cali pipe
    [~, dataflow, prmflow] = CRISrecon(calixml);
    
     % record the .corr files name
    datafile_nl(ii).output = prmflow.output;
end


