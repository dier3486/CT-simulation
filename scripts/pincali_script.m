% pin cali, after air cali

% prepare data
pindata = 'F:\data-Dier.Z\PX\bay6\200619\DetectorPosi\Temp\Pin_120KV_200MA_16x1.2_1.0_SMALLFOCAL_QFS_AIRBOWTIE_Cu.pd';

% cali xml baseline
calixmlfile = 'E:\matlab\CT\SINO\PX\Pincali_configure.xml';
calixml = readcfgfile(calixmlfile);

% output path
calioutputpath = 'F:\data-Dier.Z\PX\bay6\200619';
% input corr path
corrpath = calioutputpath;

% calibration paramters
% bad channel (shall be a corr table)
% badchannelindex = [1680	 2544   3408	4272	5136	7728	8592	9456];
badchannelindex = [];

Detectorcali = struct();
Detectorcali.pinfit_ini = [0.4      pi/2     0       0       0         0     0       0       0       0];
Detectorcali.fitselect =  [1        1        1       1       1         1     1       1       1       1];

calixml.recon.rawdata = pindata;
calixml.recon.outputpath = calioutputpath;
calixml.recon.corrpath = corrpath;
calixml.recon.pipe.Badchannel.badindex = badchannelindex;
calixml.recon.pipe.Detectorcali = Detectorcali;

fprintf('Detector calibration\n');
[~, dataflow_pin, prmflow_pin] = CRISrecon(calixml);