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
% detector_corr = [rootpath, 'system/detectorframe/detector_sample.corr'];
% detector_cfg = [rootpath, 'system/detectorframe/detector_sample.corr.xml'];
detector_corr = 'D:\matlab\ct\BCT16\detector\detector_BCT16.corr';
detector_cfg = [rootpath, 'IO\standard\detector_corr_v1.0.xml'];
detector = loadbindata(detector_corr, detector_cfg);
% clean
% detector = everything2single(detector, [], 'double');
detector.position = reshape(detector.position, [], 3);
detector.focalposition = reshape(detector.focalposition, [], 3);

% prepare
% protocal hardcode (16x0.5)
detector.position = reshape(detector.position, detector.Npixel, detector.Nslice, 3);
sliceindex = 5:20;
detector.position = reshape(detector.position(:, sliceindex, :), [], 3);
detector.Nslice = 16;

% focal, views, 
focalposition = detector.focalposition;
Nview = 1440;
% Nview = 4;
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
    [D_i, L] = intersection(focalposition, detector.position, phantom.object{iobj}, 'views-ray', viewangle, 0);
    Dmu = Dmu + D_i.*mu_i;
end

% posibility
pixel_area = 1.0;
P_air = pixel_area./(L.^2.*(pi*4));
P = exp(-Dmu).*P_air;

% P to intansity
% parameters
electric_charge = 1.602e-19;
% L0 = 1000;
% P0 = 1000;
KV = 100;
mA = 200;
mA_air = 200;
W = KV*mA;
T = 500;    % mus
gain = 0.1;

Z0 = 16384;
maxanglecode = 69120;

% Intensity
Pscale = (T*1e-6*W/electric_charge/samplekeV/1000).*gain;
Intensity = P.*Pscale + Z0;

% to put in rawdata struct
Intensity = num2cell(Intensity, 1);
readingnumber = num2cell(1:Nview, 1);
angleencoder = mod(round(viewangle./(pi*2/maxanglecode)), maxanglecode);
angleencoder = num2cell(angleencoder, 1);
rawdataversion = [1 0];
statusflag = hex2dec('8000'); 

% rawdata struct
raw(Nview) = struct();
[raw(:).Package_Version] = deal(rawdataversion);
[raw(:).Status_Flag] = deal(statusflag);
[raw(:).Reading_Number] = readingnumber{:};
[raw(:).Angle_encoder] = angleencoder{:};
[raw(:).Integration_Time] = deal(T*125);
[raw(:).KV] = deal(KV);
[raw(:).mA] = deal(mA);
[raw(:).Start_Slice] = deal(1);
[raw(:).End_Slice] = deal(detector.Nslice);
[raw(:).Raw_Data_Size] = deal(Np*3);
[raw(:).Slice_Number] = deal(detector.Nslice);
[raw(:).Raw_Data] = Intensity{:};

% rawdata output
outputpath = 'D:/data/simulation/';
% rawcfgfile = [rootpath, 'system/rawdataframe/rawdata_sample1.raw.xml'];
rawcfgfile = [rootpath, 'IO/standard/rawdata_v1.0.xml'];
outputfile = [outputpath, 'sample/rawdata_sample_project.raw'];
[raw_bin, raw_cfg] = packstruct(raw, readcfgfile(rawcfgfile), outputfile);

% offset
% skip

% air calibration
rawair = log2(P_air.*Pscale.*(mA_air/mA)) - log2(T*125);
% corr table base
aircorr_cfgfile = [rootpath, 'IO/standard/air_corr_v1.0.xml'];
aircorr_basefile = [rootpath, 'IO/standard/air_sample_v1.0.corr'];
aircorr_base = loadbindata(aircorr_basefile, aircorr_cfgfile);
% fill up 
aircorr = aircorr_base;
aircorr.ID = [0 0 1 0];
aircorr.main = rawair(:);
aircorr.mainsize = size(rawair(:),1);
aircorr.Nsection = 1;
aircorr.reference = 1.0;
aircorr.Nslice = detector.Nslice;
aircorr.Npixel = detector.Npixel;
aircorr.KV = KV;
aircorr.mA = mA_air;
% output air corr table
outaircorrfile = [outputpath, 'sample/air_sample_v1.0.corr'];
packstruct(aircorr, readcfgfile(aircorr_cfgfile), outaircorrfile);

% debug
raw2 = sparsepack(raw_bin, raw_cfg);
% read raw data
% rawdatafile = [outputpath, 'sample/rawdata_sample_project.raw'];
% raw1 = loadbindata(rawdatafile, rawcfgfile);

% recon demo (removed)

% % read air corr table
% aircorrfile = [outputpath, 'sample/air_sample_v1.0.corr'];
% aircorr1 = loadbindata(aircorrfile, aircorr_cfgfile);
% 
% % data flow
% rawhead.Angle_encoder = [raw1.Angle_encoder];
% rawhead.Reading_Number = [raw1.Reading_Number];
% rawhead.Integration_Time = [raw1.Integration_Time];
% rawhead.Reading_Number = [raw1.Reading_Number];
% rawhead.Time_Stamp = [raw1.Time_Stamp];
% rawhead.mA = [raw1.mA];
% rawhead.KV = [raw1.KV];

% log2





