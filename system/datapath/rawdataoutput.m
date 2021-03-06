function rawdataoutput(SYS, Data)
% output the rawdata and air calibration table
% [raw, aircorr]= rawdataoutput(SYS, Data);

Nw = SYS.source.Wnumber;

% values to put in rawdata 
% data version
versionstr = regexp(SYS.output.rawdataversion, '(\.\d+)|(v\d+)', 'match');
rawdataversion = [str2double(versionstr{1}(2:end)), str2double(versionstr{2}(2:end))];

% output rawdata
for iw = 1:Nw
    % put rawdata in struct
    switch rawdataversion(end-1)
        case 1
            % v1.x
            raw = rawdatastruct_v1(SYS, Data, rawdataversion, iw);
        case 2
            % v2.x
            raw = rawdatastruct_v2(SYS, Data, rawdataversion, iw);
        otherwise
            % ??
            raw = rawdatastruct_v1(SYS, Data, rawdataversion, iw);
            % pass it
    end
    
    % rawdata output
    switch SYS.output.rawdatastyle
        case {'24bit', '16bit', 'single'}
            % file name
            rawdatafile = fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} '.raw']);
            % find the format configure file
            rawcfgfile = cfgmatchrule(rawdatafile, SYS.path.IOstandard, SYS.output.rawdataversion);
            rawcfg = readcfgfile(rawcfgfile);
            % pack the data
            packstruct(raw, rawcfg, rawdatafile);
        case 'mat'
            % file name
            rawdatafile = fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} '.mat']);
            % save data
            rawdata = raw;
            save(rawdatafile, 'rawdata');
        otherwise
            % ??
            warning('Unknown style %s to save the raw data!', SYS.output.rawdatastyle);
    end
end

end

function raw = rawdatastruct_v1(SYS, Data, rawdataversion, iw)
% v1
% data cell to structure

raw = struct();
Nview = length(Data.viewangle(:));
raw(Nview) = struct();

% version ID
[raw(:).Package_Version] = deal(rawdataversion);
% status flag
[raw(:).Status_Flag] = deal(hex2dec('8000'));
% Series_Number
[raw(:).Series_Number] = deal(SYS.protocol.seriesindex);
% Shot_Number
shotnumber = num2cell(Data.shotindex, 1);
[raw(:).Shot_Number] = shotnumber{:};
% reading number
readingnumber = num2cell(1:Nview, 1);
[raw(:).Reading_Number] = readingnumber{:};
% angulation
angcode = SYS.datacollector.angulationcode;
angzero = SYS.datacollector.angulationzero;
angleencoder = mod(round(Data.viewangle./(pi*2/angcode))+angzero, angcode);
angleencoder = num2cell(angleencoder, 1);
[raw(:).Angle_encoder] = angleencoder{:};
% Integration_Time
integrationtime = round(SYS.datacollector.integrationtime*1000/SYS.datacollector.inttimeclock);
[raw(:).Integration_Time] = deal(integrationtime);
% KV
[raw(:).KV] = deal(SYS.source.KV{iw});
% mA
[raw(:).mA] = deal(SYS.source.mA{iw});
% Start_Slice
[raw(:).Start_Slice] = deal(SYS.detector.startslice);
% End_Slice
[raw(:).End_Slice] = deal(SYS.detector.endslice);
% mergescale
[raw(:).Slice_mergescale] = deal(SYS.detector.mergescale);
% Slice_Number
[raw(:).Slice_Number] = deal(max(SYS.detector.slicemerge));
% Raw_Data_Size
rawdatasize = size(Data.P{iw}, 1);
switch SYS.output.rawdatastyle
    case '16bit'
        rawdatasize = rawdatasize*2;
    case '24bit'
        rawdatasize = rawdatasize*3;
    otherwise
        if SYS.datacollector.islog2
            rawdatasize = rawdatasize*2;
        else
            rawdatasize = rawdatasize*3;
        end
end
[raw(:).Raw_Data_Size] = deal(rawdatasize);

% raw data
Raw_Data = num2cell(Data.P{iw}, 1);
[raw(:).Raw_Data] = Raw_Data{:};
    
end


function raw = rawdatastruct_v2(SYS, Data, rawdataversion, iw)
% v2
% data cell to structure

raw = struct();
Nview = length(Data.viewangle(:));
raw(Nview) = struct();
Ntube = double(SYS.source.tubenumber);

% FormatID
[raw(:).FormatID] = deal(rawdataversion);
% MachineID
% no yet
% CollectMode
log2_flag = SYS.datacollector.islog2;
rot_flag = ~strcmpi(SYS.protocol.scan, 'static');
fixangle_flag = true;
fixdose_flag = false;
CollectMode = uint16(log2_flag + rot_flag*2 + fixangle_flag*(2^2) + fixdose_flag*(2^3));
[raw(:).CollectMode] = deal(CollectMode);
% StatusFlag
[raw(:).StatusFlag] = deal(hex2dec('8000'));
% ErrorCode
[raw(:).ErrorCode] = deal(0);
% SeriesNumber
[raw(:).SeriesNumber] = deal(SYS.protocol.seriesindex);
% ShotNumber
shotnumber = num2cell(Data.shotindex, 1);
[raw(:).ShotNumber] = shotnumber{:};
% ReadingNumber
readingnumber = num2cell(1:Nview, 1);
[raw(:).ReadingNumber] = readingnumber{:};
% TimeStamp
% no yet
% AngleCode (single degree, tmp)
AngleCode = num2cell(typecast(single(Data.viewangle.*(180/pi)), 'uint32'), 1);
[raw(:).AngleCode] = AngleCode{:};
% CouchCode (single Z position, tmp)
CouchCode = num2cell(typecast(single(Data.couch(:,3)'), 'uint32'), 1);
[raw(:).CouchCode] = CouchCode{:};
% IntegrationTime (single /mus, tmp)
IntegrationTime = typecast(single(SYS.datacollector.integrationtime), 'uint32');
[raw(:).IntegrationTime] = deal(IntegrationTime);
% SetKV
[raw(:).SetKV] = deal(SYS.source.KV{iw});
% SetmA
[raw(:).SetmA] = deal(SYS.source.mA{iw});
% AECkV
[raw(:).AECkV] = deal(single(SYS.source.KV{iw}));
% AECmA
[raw(:).AECmA] = deal(single(SYS.source.mA{iw}));
% AECtime
% no yet
% StartSlice
[raw(:).StartSlice] = deal(SYS.detector.startslice);
% End_Slice
[raw(:).EndSlice] = deal(SYS.detector.endslice);
% SliceMergescale
[raw(:).SliceMergescale] = deal(SYS.detector.mergescale);
% StartPixel
if rot_flag
    TubeIndex = repmat(1:Ntube, 1, ceil(Nview/Ntube));
    TubeIndex = TubeIndex(1:Nview);
else
    TubeIndex = ones(1, Nview);
end
pixelrange = reshape(SYS.detector.pixelrange, 2, []);
StartPixel = num2cell(pixelrange(1, TubeIndex), 1); 
[raw(:).StartPixel] = StartPixel{:};
% EndPixel
EndPixel = num2cell(pixelrange(2, TubeIndex), 1); 
[raw(:).EndPixel] = EndPixel{:};
% PixelMerge
% no this function yet
[raw(:).PixelMerge] = deal(1);
% TubeIndex
TubeIndex = num2cell(TubeIndex, 1);
[raw(:).TubeIndex] = TubeIndex{:};
% CollimatorFlag
% no yet
% RawdataSize
switch SYS.output.rawdatastyle
    case '16bit'
        RawdataSize = 2;
    case '24bit'
        RawdataSize = 3;
    otherwise
        if SYS.datacollector.islog2
            RawdataSize = 2;
        else
            RawdataSize = 3;
        end
end
[raw(:).RawdataSize] = deal(RawdataSize);
% RawdataNumber
RawdataNumber = size(Data.P{iw}, 1);
[raw(:).RawdataNumber] = deal(RawdataNumber);

% rawdata
rawdata = num2cell(Data.P{iw}, 1);
[raw(:).rawdata] = rawdata{:};

end