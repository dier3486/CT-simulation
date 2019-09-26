% test script for ideal detector position
% on object
% 2D

addpath(genpath('../'));

% p = mfilename('fullpath');

% detector = load('../system/detectorframe/detectorpos_ideal_1000.mat');
detector.SID = 550;
detector.SDD = 1000;
detector.Npixel = 950;
detector.Nslice = 16;
detector.focalspot = [0, -detector.SID, 0];
detector.mid_U = 475.73;
deltafan = (54.31/180*pi)/(detector.Npixel-1);
detector.hx_ISO = detector.SID*sin(deltafan);
detector.hz_ISO = detector.hx_ISO;


righttoleft = true;
detector = idealdetector(detector, righttoleft);

detector = everything2single(detector);

% detector.slice = 1;
% detector.position = detector.position(1:950, :);
% detector.position(:, 3) = 0;

focalspot = detector.focalspot;

% object.index = 1;
% object.type = 'image2D';
% object.O = [0, 0, 0];
% object.vector = eye(3);
% object.volume = 1;

object.index = 1;
object.type = 'cylinder';
object.O = [150, 0, 0];
object.vector = [50, 0, 0; 0 50 0; 0 0 50]; 
object.volume = det(object.vector)*pi*2;

object.invV = inv(object.vector);
object.Cimage = [];

% imagepath = 'D:\Taiying\Data\testdata\Lung\';
% dcmfiles = ls([imagepath, '*.DCM']);
% object.image = [imagepath, dcmfiles(300,:)];
% info = dicominfo(object.image);
% object.Cimage = single(dicomread(info));
% 
% object.Cimage = object.Cimage+1000;
% object.Cimage(object.Cimage<-1000) = 0;
% object.Cimage = object.Cimage./1000;

% object.Cimage = zeros(512, 512);
% object.Cimage(380:410, 280:310) = 1;   


object = everything2single(object);

Nview = 1440;
% Nview = 3;
viewangle = linspace(0, pi*2, Nview+1);
viewangle = viewangle(1:end-1);

tic;
[D, L] = intersection(focalspot, detector.position, object, 'views-ray', viewangle, 0);
toc

% tic;
% D = projectioninimage(focalspot, detector.position, object.Cimage, viewangle);
% toc

% parallel projection
% parallelbeam.Np = detector.Npixel;
% parallelbeam.midchannel = detector.mid_U;
% parallelbeam.delta_d = detector.hx_ISO;
% parallelbeam.h = 1.0;

parallelbeam.Np = 950;
parallelbeam.midchannel = 476.25;
% parallelbeam.midchannel = 474.75;
parallelbeam.delta_d = 0.55;
parallelbeam.h = 500/512;
parallelbeam.viewangle = viewangle;
% parallelbeam.viewangle = mod(viewangle - pi/2, pi*2);

parallelbeam = everything2single(parallelbeam);

% fid = fopen('D:\matlab\data\simulation\20180828\td\obj.bin');
% Cimage = fread(fid,inf,'single=>single');
% fclose(fid);
% Cimage=reshape(Cimage, 512, 512);
% 
% % tic;
% % D1 = parallelprojinimage(parallelbeam, object.Cimage, '2D lineinsection');
% % toc
% 
% tic;
% D = parallelprojinimage(parallelbeam, Cimage, '2D linearinterp');
% toc