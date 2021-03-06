% step2.crs1 nonlinear#1.crs1, crosstalk#1

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PG\Nonlinearcali#1_crs2_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'E:\matlab\CT\SINO\PG\calibration\';
% namekey
namekey = 'none#1_iter2';

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
Watergoback.filter.name = 'hann';
Watergoback.filter.freqscale = 1.5;
Watergoback.span = 30;
% Watergoback.offfocal = 'deep';
% Watergoback.offfocal = 'weak';
Watergoback.offfocal = 'none';
Watergoback.offplot = true;
% crosstalk cali
crosstalkcali = struct();
crosstalkcali.Npixelpermod = 16;
crosstalkcali.Nmerge = 8;
crosstalkcali.corrversion = 'v1.11';


% debug
% datafile_nl = datafile_nl(11);
% I know the 11 is 120KV, head bowtie, small focal

% loop the protocols
Nprotocol = size(datafile_nl(:), 1);
for ii = 1:Nprotocol
    if isempty(datafile_nl(ii).filename)
        continue;
    end
    % set the values in cali xml
    calixml = calibase;
    % collimator, KV, bowtie, badchannelindex, outputpath
    for jj = 1:2
        calixml.recon{jj}.protocol.collimator = datafile_nl(ii).collimator;
        calixml.recon{jj}.protocol.bowtie = datafile_nl(ii).bowtie;
        calixml.recon{jj}.protocol.KV = datafile_nl(ii).KV;
        calixml.recon{jj}.protocol.focalsize = datafile_nl(ii).focalsize;
        calixml.recon{jj}.protocol.namekey = namekey;
        calixml.recon{jj}.outputpath = calioutputpath;
        % pipe
        calixml.recon{jj}.pipe.Air.corr = datafile_nl(ii).output.aircorr;
        calixml.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
%         calixml.recon{jj}.pipe.Crosstalk.corr = datafile_nl(ii).output.crosstalkcorr;
        calixml.recon{jj}.pipe.Nonlinear.corr = datafile_nl(ii).output.nonlinearcorr;
        calixml.recon{jj}.pipe.Offfocal = Offfocal;
        calixml.recon{jj}.pipe.Beamharden.corr = datafile_nl(ii).calitable.beamharden;
        calixml.recon{jj}.pipe.Watergoback = Watergoback;
    end
    % 1st, water 20cm off
    calixml.recon{1}.rawdata = datafile_nl(ii).filename.water200off;
    % 2nd, water 30cm off
    calixml.recon{2}.rawdata = datafile_nl(ii).filename.water300off;
    % cali
    calixml.recon{2}.pipe.crosstalkcali = crosstalkcali;
    
    % echo
    fprintf('Nonlinear Calibration #1 non-crs inter1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_nl(ii).collimator, datafile_nl(ii).bowtie, datafile_nl(ii).KV, datafile_nl(ii).focalsize);
    
    % run the cali pipe
    [~, dataflow, prmflow] = CTrecon(calixml);
    
     % record the .corr files name
    datafile_nl(ii).output.crosstalkcorr = prmflow.output.crosstalkcorr;
end


