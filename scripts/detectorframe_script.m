% detector corr table sample format
clear;

addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

% head 
detector.ID = uint8([0 0 0 0]);

detector.SID = 550;
detector.SDD = 1000;
detector.Npixel = 950;
detector.Nslice = 16;
detector.focalspot = [0, -detector.SID, 0];
detector.mid_U = 475.73;
deltafan = (54.31/180*pi)/(detector.Npixel-1);
detector.hx_ISO = detector.SID*sin(deltafan);
detector.hz_ISO = detector.hx_ISO;
detector.reserve = uint8(zeros(1, 28));

righttoleft = true;
detector = idealdetector(detector, righttoleft);

% main
righttoleft = true;
detector = idealdetector(detector, righttoleft);
detector = rmfield(detector, 'alpha_pixel');

% tail
% nothing

% cast type
detector.SID = single(detector.SID);
detector.SDD = single(detector.SDD);
detector.focalspot = single(detector.focalspot);
detector.mid_U = single(detector.mid_U);
detector.hx_ISO = single(detector.hx_ISO);
detector.hz_ISO = single(detector.hz_ISO);
detector.position = single(detector.position);
detector.Npixel = int32(detector.Npixel);
detector.Nslice = int32(detector.Nslice);

% pack to corr
[detector_corr, detector_cfg] = packstruct(detector);

% write corr (bin file)
filecorr = [rootpath, 'system/detectorframe/detector_sample.corr'];
fid = fopen(filecorr, 'w');
fwrite(fid, detector_corr, 'uint8');
fclose(fid);
% write xml format configure file
filecfg = [rootpath, 'system/detectorframe/detector_sample.corr.xml'];
root = struct();
root.detector = detector_cfg;
struct2xml(root, filecfg);

% debug
% try to read
fid = fopen(filecorr, 'r');
detector_corr2 = fread(fid, inf, 'uint8=>uint8');
fclose(fid);

filecfg2 = [rootpath, 'system/detectorframe/detector_sample2.corr.xml'];
detector_cfg2 = readcfgfile(filecfg2);
detector2 = sparsepack(detector_corr2, detector_cfg2);
% plz check if detector2==detector

% try to pack again
[detector_corr3, bincfg3] = packstruct(detector2, detector_cfg2);
% plz check if detector_corr3==detector_corr
