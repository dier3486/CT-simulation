% step1. beam harden #1

% prepare data

% the protocols to loop
toloop = struct();
% toloop.focalsize = {'small', 'large'};
toloop.focalsize = {'small'};
toloop.collimator = {'32x0.625'};
toloop.KV = [80 100 120 140];

% input data files
filepath = struct();
filepath.empty.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582870883044.0_AIR';
filepath.empty.namekey = 'empty';
filepath.empty.skiptag = 'collimator';
filepath.body.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582870883044.0_AIR';
filepath.body.namekey = 'body';
filepath.body.skiptag = 'collimator';
filepath.head.path = 'F:\data-Dier.Z\PG\bay3\DATA\1.1582870883044.0_AIR';
filepath.head.namekey = 'head';
filepath.head.skiptag = 'collimator';
% fileext
fileext = '.pd';
% get rawdata file names
datafile_bh = calidataprepare(toloop, filepath, fileext);

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PG\BHcali_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'E:\matlab\CT\SINO\PG\calibration\';

% calibration paramters
% bad channel (shall be a corr table)
badchannelindex = [];
% datamean
datamean = struct();
datamean.viewskip = 200;
% Beamhardencali
Beamhardencali = struct();
Beamhardencali.bhpolyorder = 3;
Beamhardencali.material = 'metalAl';

% debug
datafile_bh = datafile_bh(3);

% loop the protocols
Nprotocol = size(datafile_bh(:), 1);
for ii = 1:Nprotocol
    if isempty(datafile_bh(ii).filename)
        continue;
    end
    % body bowtie
    calixml_body = calibase;
    % set the values in cali xml
    calixml_body.recon{1}.rawdata = datafile_bh(ii).filename.empty;
    calixml_body.recon{2}.rawdata = datafile_bh(ii).filename.body;
    for jj = 1:2
        calixml_body.recon{jj}.protocol.collimator = datafile_bh(ii).collimator;
        calixml_body.recon{jj}.protocol.KV = datafile_bh(ii).KV;
        calixml_body.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml_body.recon{jj}.pipe.datamean = datamean;
    end
    calixml_body.recon{1}.protocol.bowtie = 'empty';
    calixml_body.recon{2}.protocol.bowtie = 'body';
    calixml_body.recon{2}.pipe.Beamhardencali = Beamhardencali;
    calixml_body.recon{2}.outputpath = calioutputpath;
    % echo
    fprintf('Beamharden Calibration #1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_bh(ii).collimator, 'body', datafile_bh(ii).KV, datafile_bh(ii).focalsize);
    % run the cali pipe
    [~, dataflow_body, prmflow_body] = CTrecon(calixml_body);
    % record the .corr files name
    datafile_bh(ii).output.beamhardencorr_body = prmflow_body.output.beamhardencorr;
    
    % head bowtie
    calixml_head = calibase;
    % set the values in cali xml
    calixml_head.recon{1}.rawdata = datafile_bh(ii).filename.empty;
    calixml_head.recon{2}.rawdata = datafile_bh(ii).filename.head;
    for jj = 1:2
        calixml_head.recon{jj}.protocol.collimator = datafile_bh(ii).collimator;
        calixml_head.recon{jj}.protocol.KV = datafile_bh(ii).KV;
        calixml_head.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml_head.recon{jj}.pipe.datamean = datamean;
    end
    calixml_head.recon{1}.protocol.bowtie = 'empty';
    calixml_head.recon{2}.protocol.bowtie = 'head';
    calixml_head.recon{2}.pipe.Beamhardencali = Beamhardencali;
    calixml_head.recon{2}.outputpath = calioutputpath;
    % echo
    fprintf('Beamharden Calibration #1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_bh(ii).collimator, 'head', datafile_bh(ii).KV, datafile_bh(ii).focalsize);
    % run the cali pipe
    [~, dataflow_head, prmflow_head] = CTrecon(calixml_head);
    % record the .corr files name
    datafile_bh(ii).output.beamhardencorr_head = prmflow_head.output.beamhardencorr;
end