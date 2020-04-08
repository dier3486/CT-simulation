% step3. nonlinear#2

% I know the datafile_nl has been done in step2.

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PG\Nonlinearcali#2_70KV_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'E:\matlab\CT\SINO\PG\calibration\';
% namekey
namekey = 'lack#2';

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
% nonlinear cali
nonlinearcali = struct();
nonlinearcali.weight = [2 1 1 1];
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
    for jj = 1:2
        calixml.recon{jj}.protocol.collimator = datafile_nl2(ii).collimator;
        calixml.recon{jj}.protocol.bowtie = datafile_nl2(ii).bowtie;
        calixml.recon{jj}.protocol.KV = datafile_nl2(ii).KV;
        calixml.recon{jj}.protocol.focalsize = datafile_nl2(ii).focalsize;
        calixml.recon{jj}.protocol.namekey = namekey;
        calixml.recon{jj}.outputpath = calioutputpath;
        % pipe
        calixml.recon{jj}.pipe.Air.corr = datafile_nl2(ii).output.aircorr;
        calixml.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml.recon{jj}.pipe.Offfocal = Offfocal;
        calixml.recon{jj}.pipe.Crosstalk.corr = datafile_nl2(ii).output.crosstalkcorr;
        calixml.recon{jj}.pipe.Beamharden.corr = datafile_nl2(ii).calitable.beamharden;
        calixml.recon{jj}.pipe.Nonlinear.corr = datafile_nl2(ii).output.nonlinearcorr;
        calixml.recon{jj}.pipe.Watergoback = Watergoback;
        
        % delete (debug)
%         calixml.recon{jj}.pipe = rmfield(calixml.recon{jj}.pipe, {'Crosstalk', 'Nonlinear'});
    end
    % 1st, water 20cm ISO
    calixml.recon{1}.rawdata = datafile_nl2(ii).filename.water200c;
    % 2nd, water 20cm off
    calixml.recon{2}.rawdata = datafile_nl2(ii).filename.water200off;
%     % 3rd, water 30cm ISO
%     calixml.recon{3}.rawdata = datafile_nl2(ii).filename.water300c;
%     % 4th, water 30cm off
%     calixml.recon{4}.rawdata = datafile_nl2(ii).filename.water300off;
    % nonlinear cali
    calixml.recon{2}.pipe.nonlinearcali = nonlinearcali;

    % echo
    fprintf('Nonlinear Calibration #2 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_nl2(ii).collimator, datafile_nl2(ii).bowtie, datafile_nl2(ii).KV, datafile_nl2(ii).focalsize);
    
    % run the cali pipe
    [~, dataflow, prmflow] = CTrecon(calixml);
    
     % record the .corr files name
    datafile_nl2(ii).output.nonlinearcorr = prmflow.output.nonlinearcorr;
end


