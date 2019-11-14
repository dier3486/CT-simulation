% corr table basic format
clear;

addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

% air corr
air.ID = uint8([0 0 0 0]);
air.Npixel = int32(950);
air.Nslice = int32(16);
air.startslice = int32(5);
air.endslice = int32(20);
air.slicemerge = int32(1);
air.focalspot = int32(1);
air.KV = int32(120);
air.mA = int32(30);
air.bowtie = int32(0);
air.rotationspeed = single(1.0);
air.Nsection = int32(20);
air.firstangle = single(0);
air.mainsize = int32(air.Npixel*air.Nslice*air.Nsection);
air.reserve = uint8([]);
air.reference = single(zeros(air.Nsection, 1));
air.main = single(zeros(air.Npixel*air.Nslice*air.Nsection, 1));
% pack to corr
[air_corr1, air_cfg1] = packstruct(air);

% fix format
air_cfg = air_cfg1;
air_cfg.size = [];
air_cfg.reserve.number = 512-air_cfg1.reference.offset;
air_cfg.reference.offset = 512;
% air_cfg.reference.number = double(air.Nsection);
air_cfg.reference.number = '$.Nsection';
air_cfg.main.offset = 1024;
air_cfg.main.number = '$.mainsize';

% write xml format configure file
filecfg = [rootpath, 'IO/standard/air_corr_v1.0.xml'];
root = struct();
root.air = air_cfg;
struct2xml(root, filecfg);

% repack to corr
% write corr (bin file)
filecorr = [rootpath, 'IO/standard/air_sample_v1.0.corr'];
[air_corr, air_cfg2] = packstruct(air, air_cfg, filecorr);


