% detector corr table sample format
addpath(genpath('../'));
rootpath = 'D:/matlab/CTsimulation/';

raw = struct();
% head 
raw.start = uint32(0);
raw.Data_Flag = uint8(1);
raw.Sample_Mode = uint8(0);
raw.Package_Version = uint8([0 0]);
raw.Status_Flag = uint16(0);
raw.Error_Code = uint16(0);
raw.Series_Number = uint16(0);
raw.Shot_Number = uint16(0);
raw.Reading_Number = uint32(0);
raw.Time_Stamp = uint32(0);
raw.Angle_encoder = uint32(0);
raw.Table_encoder = uint32(0);
raw.Integration_Time = uint32(0);
raw.reserved1 = uint8([0 0 0]);
raw.Focal_Spot = uint8(0);
raw.KV = uint16(0);
raw.mA = uint16(0);
raw.Filament = uint32(0);
raw.Start_Slice = uint8(1);
raw.End_Slice = uint8(16);
raw.Slice_merge = uint8(1);
raw.Slice_Number = uint8(16);
raw.Raw_Data_Size = uint32(950*16*3);
raw.Raw_Data_Address = uint32(0);
raw.reserved2 = uint8([0 0 0 0]);

% main
raw.Raw_Data = uint8(1:950*16*3);

Nview = 100;
raw(Nview) = raw(1);
[raw(:)] = raw(1);

% tail
% nothing

% pack to bin
[raw_bin, raw_cfg] = packstruct(raw);

% write raw (bin file)
fileraw = [rootpath, 'system/rawdataframe/rawdata_sample.raw'];
fid = fopen(fileraw, 'w');
fwrite(fid, raw_bin, 'uint8');
fclose(fid);
% write xml format configure file
filecfg = [rootpath, 'system/rawdataframe/rawdata_sample0.raw.xml'];
root = struct();
root.raw = raw_cfg;
struct2xml(root, filecfg);

% debug
% try to read
fid = fopen(fileraw, 'r');
raw_bin2 = fread(fid, inf, 'uint8=>uint8');
fclose(fid);

filecfg2 = [rootpath, 'system/rawdataframe/rawdata_sample1.raw.xml'];
raw_cfg2 = readcfgfile(filecfg2);
raw2 = sparsepack(raw_bin2, raw_cfg2);

raw_bin3 = packstruct(raw2, raw_cfg2);

