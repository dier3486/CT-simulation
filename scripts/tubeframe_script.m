% detector corr table sample format
clear;

addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

% head 
tube.ID = uint8([0 0 0 0]);
tube.focalsize = single([0.9 0.7 1.4 1.4]);
tube.anodeangle = single(7.0);
tube.KVnumber = uint32(4);
tube.KVtag = single([80 100 120 140]);
tube.Nsample = uint32(140);

% tube spectrum
% load data
datapath = 'D:\matlab\CTsimulation\physics\tube\';
spect80 = csvread([datapath, '80kVp.csv'], 1, 0);
spect100 = csvread([datapath, '100kVp.csv'], 1, 0);
spect120 = csvread([datapath, '120kVp.csv'], 1, 0);
spect140 = csvread([datapath, '140kVp.csv'], 1, 0);
% I know they are in same size \o/
tube.main = single([spect80 spect100 spect120 spect140]);

% pack to corr
[tube_corr, tube_cfg] = packstruct(tube);

% fix cfg
tube_cfg.size = [];
tube_cfg.main.offset = 512;
tube_cfg.main.number = '$.Nsample*$.KVnumber*2';

% write xml format configure file
filecfg = [rootpath, 'IO/standard/tube_corr_v1.0.xml'];
root = struct();
root.tube = tube_cfg;
struct2xml(root, filecfg);

% write bin file
binfile = [datapath, 'tube_spectrumdata_v1.0.corr'];
[tube_corr2, tube_cfg2] = packstruct(tube, tube_cfg, binfile);

