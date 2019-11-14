% detector corr table sample format
clear;

addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

% head 
bowtie.ID = uint8([0 0 0 0]);
bowtie.box = single([170, 20, 55]);
bowtie.focaltobottom = single(159.2);
bowtie.Ncurve = uint32(2);
bowtie.Nsample = uint32(0);

% bowtie curve
% load data
datafile = 'D:\matlab\ct\BCT16\collimation\bowtie full 170mm.csv';
bowtie.main = single(csvread(datafile));

% pack to corr
[bowtie_corr, bowtie_cfg] = packstruct(bowtie);

% fix cfg
bowtie.Nsample = uint32(size(bowtie.main,1));
bowtie_cfg.size = [];
bowtie_cfg.main.offset = 512;
bowtie_cfg.main.number = '$.Nsample*($.Ncurve+1)';

% write xml format configure file
filecfg = [rootpath, 'IO/standard/bowtie_corr_v1.0.xml'];
root = struct();
root.bowtie = bowtie_cfg;
struct2xml(root, filecfg);

% write bin file
binfile = 'D:\matlab\ct\BCT16\collimation\bowtie_geometry_v1.0.corr';
[bowtie_corr2, bowtie_cfg2] = packstruct(bowtie, bowtie_cfg, binfile);

