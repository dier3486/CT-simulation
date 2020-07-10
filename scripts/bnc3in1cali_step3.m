% step3. nonlinear#2

% I know the datafile_nl has been done in step2.

% cali xml baseline
% calixmlfile = 'E:\matlab\CT\SINO\TM\Nonlinearcali#2_configure.xml';
calixmlfile = 'E:\matlab\CT\SINO\PG\Nonlinearcali#2_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'F:\data-Dier.Z\PG\bay3\20200703\';
% namekey
% namekey = 'none#2';
% input corr path (to looking for .corr files in this folder)
inputcorrpath = calioutputpath;

% calibration paramters
% bad channel (shall be a corr table)
badchannelindex = [];
% badchannelindex = [1680	2544	3408	4272	5136	7728	8592	9456];
% off-focal corr (shall be a corr table)
Offfocal = struct();
% Offfocal.offintensity = [0.0007 0.000];
% Offfocal.offwidth = [110 0];
% Offfocal.offedge = [0.6 0.6];
% Offfocal.ratescale = [1.0 0.8];
Offfocal.offintensity = [0.005 0.001];
Offfocal.offwidth = [65 95];
Offfocal.offedge = [0.6 0.6];
Offfocal.ratescale = [0.8 0.8];
Offfocal.crossrate = 0.75;
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
nonlinearcali.weight = [3 1 1 1];
nonlinearcali.minptrans = 16;
nonlinearcali.corrversion = 'v1.11';



datafile_nl2 = datafile_nl;
% % debug
% datafile_nl = datafile_nl(11);
% I know the 11 is 120KV, head bowtie, small focal

% loop the protocols
Nprotocol = size(datafile_nl2(:), 1);
for ii = 1:Nprotocol
    if isempty(datafile_nl2(ii).filename)
        continue;
    end
    % set the values in cali xml
    calixml = calibase;
    % collimator, KV, bowtie, badchannelindex, outputpath, and pipe line
    for jj = 1:4
        calixml.recon{jj}.outputpath = calioutputpath;
        calixml.recon{jj}.corrpath = inputcorrpath;
        calixml.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml.recon{jj}.pipe.Offfocal = structmerge(Offfocal, calixml.recon{jj}.pipe.Offfocal);
        calixml.recon{jj}.pipe.Idealwater = structmerge(Watergoback, calixml.recon{jj}.pipe.Idealwater);
    end
    % 1st, water 20cm ISO
    calixml.recon{1}.rawdata = datafile_nl2(ii).filename.water200c;
    calixml.recon{1}.protocol.namekey = 'SMALLWATER_CENTER';
    % 2nd, water 20cm off
    calixml.recon{2}.rawdata = datafile_nl2(ii).filename.water200off;
    calixml.recon{2}.protocol.namekey = 'SMALLWATER_OFFCEN';
    % 3rd, water 30cm ISO
    calixml.recon{3}.rawdata = datafile_nl2(ii).filename.water300c;
    calixml.recon{3}.protocol.namekey = 'LARGEWATER_CENTER';
    % 4th, water 30cm off
    calixml.recon{4}.rawdata = datafile_nl2(ii).filename.water300off;
    calixml.recon{4}.protocol.namekey = 'LARGEWATER_OFFCEN';
    % nonlinear cali
    calixml.recon{4}.pipe.nonlinearcali = structmerge(nonlinearcali, calixml.recon{4}.pipe.nonlinearcali);

    % echo
    fprintf('Nonlinear Calibration #2 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_nl2(ii).collimator, datafile_nl2(ii).bowtie, datafile_nl2(ii).KV, datafile_nl2(ii).focalsize);
    
    % run the cali pipe
    [~, dataflow, prmflow] = CRISrecon(calixml);
    
     % record the .corr files name
    datafile_nl2(ii).output.nonlinearcorr = prmflow.output.nonlinearcorr;
end


