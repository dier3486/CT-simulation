% a projection script for single energy
clear;
addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

% single energy
samplekeV = 65;

% configure phantom
phantom_cfg = 'sample_phantom.json';
phantom = phantomconfigure(phantom_cfg);
% configure phantom's material
phantom = materialconfigure(phantom, samplekeV);

% load detector
detector_corr = [rootpath, 'system/detectorframe/detector_sample.corr'];
detector_cfg = [rootpath, 'system/detectorframe/detector_sample.corr.xml'];
detector = loadbindata(detector_corr, detector_cfg);
% clean
% detector = everything2single(detector, [], 'double');
detector.position = reshape(detector.position, [], 3);

% prepare
focalspot = detector.focalspot;
Nview = 1440;
% Nview = 3;
viewangle = linspace(0, pi*2, Nview+1);
viewangle = viewangle(1:end-1);
Np = detector.Npixel * detector.Nslice;

% projection
Dmu = zeros(Np, Nview);
for iobj = 1:phantom.Nobject
    mu_i = phantom.object{iobj}.material.mu_total;
    parentobj = phantom.object_tree(iobj);
    if parentobj>0
        mu_i = mu_i - phantom.object{parentobj}.material.mu_total;
    end
    [D_i, L] = intersection(focalspot, detector.position, phantom.object{iobj}, 'views-ray', viewangle, 0);
    Dmu = Dmu + D_i.*mu_i;
end

% posibility
pixel_area = 1.0;
P = exp(-Dmu).*pixel_area./(L.^2.*(pi*4));

% P to intansity
% parameters
electric_charge = 1.602e-19;
% L0 = 1000;
% P0 = 1000;
kV = 100;
mA = 200;
W = kV*mA;
T = 500;    % mus
gain = 0.1;

Z0 = 16384;
maxanglecode = 69120;

Intensity = P.*(T*1e-6*W/electric_charge/samplekeV/1000).*gain + Z0;
Intensity = num2cell(Intensity, 1);
readingnumber = num2cell(1:Nview, 1);
angleencoder = mod(round(viewangle./(pi*2/maxanglecode)), maxanglecode);
angleencoder = num2cell(angleencoder, 1);

% rawdata
raw(Nview) = struct();
[raw(:).Reading_Number] = readingnumber{:};
[raw(:).Angle_encoder] = angleencoder{:};
[raw(:).Integration_Time] = deal(T*125);
[raw(:).KV] = deal(kV);
[raw(:).mA] = deal(mA);
[raw(:).Start_Slice] = deal(1);
[raw(:).End_Slice] = deal(detector.Nslice);
[raw(:).Raw_Data_Size] = deal(Np*3);
[raw(:).Slice_Number] = deal(detector.Nslice);
[raw(:).Raw_Data] = Intensity{:};

% rawdata format cfg
cfgfile = [rootpath, 'system/rawdataframe/rawdata_sample1.raw.xml'];
raw_cfg = readcfgfile(cfgfile);

% pack
raw_bin = packstruct(raw, raw_cfg);

% output
outputpath = 'D:/data/simulation/';
outputfile = [outputpath, 'sample/rawdata_sample_project.raw'];
fid = fopen(outputfile, 'w');
fwrite(fid, raw_bin, 'uint8');
fclose(fid);

% debug
raw2 = sparsepack(raw_bin, raw_cfg);



