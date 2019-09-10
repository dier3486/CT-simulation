function SYS = systemcfgsample(detector, phantom)
% return a default CT system config, just a sample

if nargin<1
    detector = idealdetector();
end

if nargin<2
    phantom = phantomcfgsample();
end

% path
SYS.path.main = 'D:\matlab\CTsimulation\';
SYS.path.matter = 'D:\matlab\CTsimulation\physics\matter\';
SYS.path.output = 'D:\matlab\data\simulation\';

% world
SYS.world.samplekeV_range = 1:150;
SYS.world.refrencekeV = 60;

% detector
SYS.detector = detector;
SYS.detector.reponse = 1.0;
% ASG
SYS.detector.ASG = [];
% fliter on detector
SYS.detector.filter = [];

% source
SYS.source.focalposition = [0, -detector.SID, 0];
SYS.source.spectrum = []; % data path

% collimation
SYS.collimation.bowtie = []; % data path
SYS.collimation.filter = [];
SYS.collimation.blades = [];

% protocal
SYS.protocal{1}.method = 'geometry';
SYS.protocal{1}.scan = 'Axial';
SYS.protocal{1}.collimator = '16x0.55';
SYS.protocal{1}.correction = {};

% scatter
SYS.scatter = [];

% phantom
SYS.phantom = phantom;
