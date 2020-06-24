% step1. beam harden #1

% prepare data

% the protocols to loop
toloop = struct();
% toloop.focalsize = {'small', 'large'};
toloop.focalsize = {'SMALL'};
toloop.focaltype = {'QFS'};
toloop.collimator = {'16x1.2'};
toloop.KV = [120];

bhdatapath = 'F:\data-Dier.Z\PX\bay6\BH\';
% input data files
filepath = struct();
filepath.empty.path = bhdatapath;
filepath.empty.namekey = {'HardenA', 'AIRBOWTIE'};
% filepath.empty.skiptag = 'collimator';
filepath.body.path = bhdatapath;
filepath.body.namekey = {'HardenA', 'LARGEBOWTIE'};
% filepath.body.skiptag = 'collimator';
filepath.head.path = bhdatapath;
filepath.head.namekey = {'HardenA', 'SMALLBOWTIE'};
% filepath.head.skiptag = 'collimator';
% fileext
fileext = '.pd';
% get rawdata file names
datafile_bh = calidataprepare(toloop, filepath, fileext);

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PX\BHcali_configure.xml';
calibase = readcfgfile(calixmlfile);

% output path
calioutputpath = 'F:\data-Dier.Z\PX\bay6\';

% calibration paramters
% bad channel (shall be a corr table)
% badchannelindex = [1680	2544	3408	4272	5136	7728	8592	9456];
badchannelindex = [];
% datamean
datamean = struct();
datamean.viewskip = 200;
% Beamhardencali
Beamhardencali = struct();
Beamhardencali.bhpolyorder = 3;
Beamhardencali.material = 'teflon';
Beamhardencali.corrversion = 'v1.11';
Beamhardencali.toplot = true;

% debug
% datafile_bh = datafile_bh(3);

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
%         calixml_body.recon{jj}.protocol.collimator = datafile_bh(ii).collimator;
%         calixml_body.recon{jj}.protocol.KV = datafile_bh(ii).KV;
%         calixml_body.recon{jj}.protocol.focalsize = datafile_bh(ii).focalsize;
        calixml_body.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml_body.recon{jj}.pipe.datamean = datamean;
        calixml_body.recon{jj}.outputpath = calioutputpath;
    end
%     calixml_body.recon{1}.protocol.bowtie = 'empty';
%     calixml_body.recon{2}.protocol.bowtie = 'body';
    calixml_body.recon{2}.pipe.Beamhardencali = Beamhardencali;
    
    % echo
    fprintf('Beamharden Calibration #1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_bh(ii).collimator, 'body', datafile_bh(ii).KV, datafile_bh(ii).focalsize);
    % run the cali pipe
    [~, dataflow_body, prmflow_body] = CRISrecon(calixml_body);
    % record the .corr files name
    datafile_bh(ii).output.beamhardencorr_body = prmflow_body.output.beamhardencorr;
    
    % head bowtie
    calixml_head = calibase;
    % set the values in cali xml
    calixml_head.recon{1}.rawdata = datafile_bh(ii).filename.empty;
    calixml_head.recon{2}.rawdata = datafile_bh(ii).filename.head;
    for jj = 1:2
%         calixml_head.recon{jj}.protocol.collimator = datafile_bh(ii).collimator;
%         calixml_head.recon{jj}.protocol.KV = datafile_bh(ii).KV;
%         calixml_head.recon{jj}.protocol.focalsize = datafile_bh(ii).focalsize;
        calixml_head.recon{jj}.pipe.Badchannel.badindex = badchannelindex;
        calixml_head.recon{jj}.pipe.datamean = datamean;
        calixml_head.recon{jj}.outputpath = calioutputpath;
    end
%     calixml_head.recon{1}.protocol.bowtie = 'empty';
%     calixml_head.recon{2}.protocol.bowtie = 'head';
    calixml_head.recon{2}.pipe.Beamhardencali = Beamhardencali;
    % echo
    fprintf('Beamharden Calibration #1 for: %s, %s Bowtie, %d KV, %s Focal\n', ...
        datafile_bh(ii).collimator, 'head', datafile_bh(ii).KV, datafile_bh(ii).focalsize);
    % run the cali pipe
    [~, dataflow_head, prmflow_head] = CRISrecon(calixml_head);
    % record the .corr files name
    datafile_bh(ii).output.beamhardencorr_head = prmflow_head.output.beamhardencorr;
end