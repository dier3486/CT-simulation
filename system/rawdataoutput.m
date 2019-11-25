function [raw, aircorr]= rawdataoutput(SYS, Data)
% output the raw{iw} data

Nw = SYS.source.Wnumber;
raw = cell(1, Nw);
aircorr = cell(1, Nw);

% values to put in rawdata and/or air_corr
% data version
versionstr = regexp(SYS.output.rawdataversion, '(\.\d+)|(v\d+)', 'match');
rawdataversion = [str2double(versionstr{1}(2:end)), str2double(versionstr{2}(2:end))];
% status
statusflag = hex2dec('8000');
% Series_Number
seriesnumber = SYS.protocol.series_index;
% Shot_Number
shotnumber = num2cell(Data.shotindex, 1);
% Nview
Nview = length(Data.viewangle(:));
% reading number
readingnumber = num2cell(1:Nview, 1);
% angulation
angcode = SYS.datacollector.angulationcode;
angzero = SYS.datacollector.angulationzero;
angleencoder = mod(round((Data.viewangle-pi/2)./(pi*2/angcode))+angzero, angcode);
angleencoder = num2cell(angleencoder, 1);
% Integration_Time
integrationtime = round(SYS.datacollector.integrationtime*1000/SYS.datacollector.inttimeclock);
% Start_Slice
startslice = SYS.detector.startslice;
% End_Slice
endslice = SYS.detector.endslice;
% Slice_merge
slicemerge = SYS.detector.mergescale;
% Slice_Number
slicenumber = max(SYS.detector.slicemerge);
% Raw_Data_Size
rawdatasize = SYS.detector.Npixel*slicenumber*3;    % 24bit
% Npixel
Npixel = SYS.detector.Npixel;

% output rawdata
for iw = 1:Nw
    % raw data to put in struct
    raw_Data = num2cell(Data.P{iw}, 1);
    % values to put in struct depending on KV mA
    KV = SYS.source.KV{iw};
    mA = SYS.source.mA{iw};
    % rawdata struct
    raw{iw} = struct();
    raw{iw}(Nview) = struct();
    [raw{iw}(:).Package_Version] = deal(rawdataversion);
    [raw{iw}(:).Status_Flag] = deal(statusflag);
    [raw{iw}(:).Series_Number] = deal(seriesnumber);
    [raw{iw}(:).Shot_Number] = shotnumber{:};
    [raw{iw}(:).Reading_Number] = readingnumber{:};
    [raw{iw}(:).Angle_encoder] = angleencoder{:};
    [raw{iw}(:).Integration_Time] = deal(integrationtime);
    [raw{iw}(:).KV] = deal(KV);
    [raw{iw}(:).mA] = deal(mA);
    [raw{iw}(:).Start_Slice] = deal(startslice);
    [raw{iw}(:).End_Slice] = deal(endslice);
    [raw{iw}(:).Slice_merge] = deal(slicemerge);
    [raw{iw}(:).Slice_Number] = deal(slicenumber);
    [raw{iw}(:).Raw_Data_Size] = deal(rawdatasize);
    [raw{iw}(:).raw_Data] = raw_Data{:};
    
    % rawdata output
    % file name
    rawdatafile = [SYS.output.path SYS.output.files.rawdata{iw} '_' SYS.output.rawdataversion '.raw'];
    % find the format configure file
    rawcfgfile = cfgmatchrule(rawdatafile, SYS.path.IOstandard, SYS.output.rawdataversion);
    rawcfg = readcfgfile(rawcfgfile);
    % pack the data
    packstruct(raw{iw}, rawcfg, rawdatafile);
end

% output air corr
% corr table baseline
aircorr_basefile = [SYS.path.IOstandard, 'air_sample_v1.0.corr'];
if exist(aircorr_basefile, 'file')
    aircorr_base = loaddata(aircorr_basefile, SYS.path.IOstandard);
else
    % empty baseline
    aircorr_base = struct();
end
for iw = 1:Nw
    % copy from base
    aircorr{iw} = aircorr_base;
    % values to put in struct depending on KV mA
    KV = SYS.source.KV{iw};
    mA_air = SYS.source.mA_air{iw};
    % air main
    aircorr{iw}.ID = [0 0 1 0];
    aircorr{iw}.main = single(Data.Pair{iw}(:));
    aircorr{iw}.mainsize = length(aircorr{iw}.main);
    aircorr{iw}.Nsection = 1;
    aircorr{iw}.reference = 1.0;
    aircorr{iw}.Nslice = slicenumber;
    aircorr{iw}.Npixel = Npixel;
    aircorr{iw}.KV = KV;
    aircorr{iw}.mA = mA_air;
    
    % output air corr table
	aircorrfile = [SYS.output.path SYS.output.files.aircorr{iw} '_v1.0' '.corr'];
    aircfgfile = cfgmatchrule(aircorrfile, SYS.path.IOstandard, 'v1.0');
    aircfg = readcfgfile(aircfgfile);
    packstruct(aircorr{iw}, aircfg, aircorrfile);
end

end