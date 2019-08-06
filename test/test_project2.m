% test script for ideal detector position
% on grid image
% 2D

addpath(genpath('../'));

detector = load('../system/detectorframe/detectorpos_ideal_1000.mat');

focalspot = [0, -detector.SSD, 0];

object.index = 1;
object.type = 'image3D';
object.O = [0, 0, 0];
object.vector = eye(3);
object.volume = 1;

object.invV = inv(object.vector);

imagepath = 'D:\Taiying\Data\testdata\Lung\';
dcmfiles = ls([imagepath, '*.DCM']);
object.image = [imagepath, dcmfiles(300,:)];
info = dicominfo(object.image);
object.Cimage = single(dicomread(info));

object.Cimage = object.Cimage+1000;
object.Cimage(object.Cimage<-1000) = 0;
object.Cimage = object.Cimage./1000;

% Nview = 1024;
Nview = 3;
viewangle = linspace(0, pi, Nview+1);
viewangle = viewangle(1:end-1);

% tic;
% [D, L] = intersection(focalspot, detector.position, object, 'views-ray', viewangle, 0);
% toc

tic;
D = projectioninimage(focalspot, detector.position, object.Cimage, viewangle);
toc
