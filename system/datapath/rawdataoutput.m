function rawdataoutput(SYS, Data)
% output the rawdata and air calibration table
% [raw, aircorr]= rawdataoutput(SYS, Data);

Nw = SYS.source.Wnumber;
iblock = Data.iblock;

% values to put in rawdata 
% data version
versionstr = regexp(SYS.output.rawdataversion, '(\.\d+)|(v\d+)', 'match');
rawdataversion = [str2double(versionstr{1}(2:end)), str2double(versionstr{2}(2:end))];

% output rawdata
for iw = 1:Nw
    % put rawdata in struct
    switch rawdataversion(end-1)
        case {0, 1}
            % v1.x or 0.x
            raw = simurawdatastruct(SYS, Data, rawdataversion, iw);
        case 2
            % v2.x
            raw = simurawdatastructV2(SYS, Data, rawdataversion, iw);
        otherwise
            % ??
            raw = simurawdatastruct(SYS, Data, rawdataversion, iw);
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
            if iblock == 1
                packstruct(raw, rawcfg, rawdatafile);
            else
                packstruct(raw, rawcfg, rawdatafile, 'a');
            end
        case 'mat'
            % file name
            if iblock == 1
                rawdatafile = fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} '.mat']);
            else
                exname = sprintf('.m%02d', iblock-1);
                rawdatafile = fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} exname]);
            end
            % save data
            rawdata = raw;
            save(rawdatafile, 'rawdata');
        otherwise
            % ??
            warning('Unknown style %s to save the raw data!', SYS.output.rawdatastyle);
    end
end

end

function raw = simurawdatastruct(SYS, Data, rawdataversion, iw)
% v1
% data cell to structure

raw = struct();
Nview = length(Data.viewangle(:));
startreading = Data.startreading;
raw(Nview) = struct();

% version ID
[raw(:).PackageVersion] = deal(rawdataversion);
% status flag
[raw(:).StatusFlag] = deal(hex2dec('8000'));
% SeriesNumber
[raw(:).SeriesNumber] = deal(SYS.protocol.seriesindex);
% ShotNumber
shotnumber = num2cell(Data.shotindex, 1);
[raw(:).ShotNumber] = shotnumber{:};
% reading number
% readingnumber = num2cell(1:Nview, 1);
readingnumber = num2cell((0:Nview-1)+startreading, 1);
[raw(:).ReadingNumber] = readingnumber{:};
% angulation
angcode = SYS.datacollector.angulationcode;
angzero = SYS.datacollector.angulationzero;
angleencoder = mod(round(Data.viewangle./(pi*2/angcode))+angzero, angcode);
angleencoder = num2cell(angleencoder, 1);
[raw(:).AngleEncoder] = angleencoder{:};
% IntegrationTime
integrationtime = round(SYS.datacollector.integrationtime*1000/SYS.datacollector.inttimeclock);
[raw(:).IntegrationTime] = deal(integrationtime);
% TimeStamp
TimeStamp = zeros(1, Nview);
stamp0 = floor(rand(1)*2^32);
stampperview = SYS.protocol.rotationspeed*1e9 / SYS.protocol.viewperrot / SYS.datacollector.inttimeclock;
stamp_rot0 = stamp0;
for ii = 1 : ceil(SYS.protocol.rotationnumber)
    N0 = (ii-1)*SYS.protocol.viewperrot;
    Nv = min(SYS.protocol.viewperrot, Nview - N0);
    TimeStamp(N0+(1:Nv)) = mod((1:Nv).*stampperview + stamp_rot0, 2^32);
    stamp_rot0 = TimeStamp(N0+Nv);
end
TimeStamp = num2cell(uint32(TimeStamp), 1);
[raw(:).TimeStamp] = TimeStamp{:};
% TableEncoder
movercode = mod(Data.couch(:,3)'./SYS.datacollector.moverlength.*2^SYS.datacollector.movercode, 2^SYS.datacollector.movercode);
movercode = round(movercode.*SYS.datacollector.moveruppersample);
movercode = num2cell(uint32(movercode), 1);
[raw(:).TableEncoder] = movercode{:};
% TableGear
if isfield(Data, 'couchgear')
    couchgear = num2cell(int8(Data.couchgear), 1);
    [raw(:).TableGear] = couchgear{:};
end

% KV
[raw(:).KV] = deal(SYS.source.KV{iw});
% mA
[raw(:).mA] = deal(SYS.source.mA{iw});
% StartSlice
[raw(:).StartSlice] = deal(SYS.detector.startslice);
% EndSlice
[raw(:).EndSlice] = deal(SYS.detector.endslice);
% mergescale
[raw(:).SliceMergescale] = deal(SYS.detector.mergescale);
% SliceNumber
[raw(:).SliceNumber] = deal(max(SYS.detector.slicemerge));
% RawDataSize
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
[raw(:).RawDataSize] = deal(rawdatasize);

% raw data
RawData = num2cell(Data.P{iw}, 1);
[raw(:).RawData] = RawData{:};
    
end


function raw = rawdatastruct_v2(SYS, Data, rawdataversion, iw)
% v2
% data cell to structure

raw = struct();
Nview = length(Data.viewangle(:));
startreading = Data.startreading;
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
readingnumber = num2cell((0:Nview-1)+startreading, 1);
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
% EndSlice
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