function SYS_cfg = systemcfgsample(rootpath)
% return a default CT system config, just a sample

if nargin<1
    rootpath = '';
end
if isempty(rootpath)
    rootpath = pwd;
end

% path
SYS_cfg.path.main = rootpath;
SYS_cfg.path.matter = '~\physics\matter\';
SYS_cfg.path.IOstandard = '~\IO\standard\';
SYS_cfg.path.systemdata = '~\system\';

% world
SYS_cfg.world.elementsdata = '$matter\elements\';
SYS_cfg.world.materialdata = '$matter\material\';
SYS_cfg.world.samplekeV_range = [5, 150];
SYS_cfg.world.samplekeV_step = 0.5;
SYS_cfg.world.referencekeV = 60;
SYS_cfg.world.water.material = 'water';

% detector
SYS_cfg.detector.frame_base = '$systemdata\detectorframe\detector_sample_v1.0.corr';
SYS_cfg.detector.frame_extra = [];
SYS_cfg.detector.reponse = 1.0;
% ASG (on detector)
SYS_cfg.detector.ASG = [];
% fliter (on detector)
SYS_cfg.detector.filter = [];

% source
SYS_cfg.source.focalposition = [0 -568 0];
SYS_cfg.source.focaldistort = 0;
% SYS_cfg.source.focalsize = [0.7, 1.0];
SYS_cfg.source.tubedata = '~\physics\tube\tube_spectrumdata_v1.0.corr';

% collimation
SYS_cfg.collimation.bowtie.bowtiedata = '$systemdata\bowtieframe\bowtie_sample_v1.0.corr';
SYS_cfg.collimation.bowtie.material = 'Teflon';
SYS_cfg.collimation.filter{1}.thickness = 2.0;
SYS_cfg.collimation.filter{1}.material = 'metalAl';
SYS_cfg.collimation.filter{2}.thickness = 1.0;
SYS_cfg.collimation.filter{2}.material = 'metalTi';
SYS_cfg.collimation.blades.blasesdata = '';

% DCB
SYS_cfg.datacollector.angulationcode = 69120;
SYS_cfg.datacollector.angulationzero = 0;
SYS_cfg.datacollector.DBBzero = 16384;      %Z0
SYS_cfg.datacollector.DBBgain = 0.1;
SYS_cfg.datacollector.inttimeclock = 8; % ns

% console
SYS_cfg.console.protocoltrans = '';
SYS_cfg.console.dicomdictionary = '';

% simulation method
SYS_cfg.simulation.project = 'Geometry';
SYS_cfg.simulation.spectrum = 'Single';  % or Continue
SYS_cfg.simulation.detectsample = 1;
SYS_cfg.simulation.focalsample = 1;
SYS_cfg.simulation.quantumnoise = 0;
SYS_cfg.simulation.offfocal = 0;
SYS_cfg.simulation.scatter = 0;

% scatter paramters (TBC)
SYS_cfg.scatter = [];

% output
SYS_cfg.output.path = 'D:\matlab\data\simulation\';
SYS_cfg.output.namekey = '';
SYS_cfg.output.namerule = 'simple';
SYS_cfg.output.rawdataversion = 'v1.0';
SYS_cfg.output.corrtable = 'air, beamharden';
SYS_cfg.output.corrversion = [];    % default


