% pin cali, after air cali

% prepare data
pindata = 'F:\data-Dier.Z\PG\bay4\20200612\pin';

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PG\Pincali_configure.xml';
calixml = readcfgfile(calixmlfile);

% output path
calioutputpath = 'F:\data-Dier.Z\PG\bay4\20200612\';
% input corr path
corrpath = calioutputpath;

% calibration paramters
% bad channel (shall be a corr table)
badchannelindex = [1680	 2544   3408	4272	5136	7728	8592	9456];
% badchannelindex = [];

calixml.recon.rawdata = pindata;
calixml.recon.outputpath = calioutputpath;
calixml.recon.corrpath = corrpath;
calixml.recon.pipe.Badchannel.badindex = badchannelindex;

fprintf('Detector calibration\n');
[~, dataflow_pin, prmflow_pin] = CRISrecon(calixml);