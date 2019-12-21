function CRIS2rawdata(pdfile, outpath, CRISpath, cfgpath, versionflag)
% trans .pd rawdata to .raw
% CRIS2rawdata(pdfile, outpath, CRISpath, cfgpath, versionflag);
% INPUT:
%   pdfile          filename for a cell of filenames
%   outpath         output path
%   CRISpath        external IO platform
%   cfgpath         ~\IO\standard\
%   versionflag     version, e.g. 'v1.0'

% default inputs
if nargin<2
    outpath = '';
end
if nargin<3
    CRISpath = '';
end
if nargin<4
    cfgpath = '';
end
if nargin<5
    versionflag = '';
end

% loop multi files
if iscell(pdfile)
    % recurse
    for ii = 1:length(pdfile)
         CRIS2rawdata(pdfile{ii}, outpath, CRISpath, cfgpath, versionflag)
    end
    return;
end

% parts in input file
[pdpath, pdname, ~] = fileparts(pdfile);

if isempty(outpath)
    % default output path is the orig path
    outpath = pdpath;
end

if isempty(CRISpath)
    % default CRIS path
    CRISpath = 'E:\matlab\CRIS';
    % I know it is
end

if isempty(versionflag)
    versionflag = 'v1.0';
end

% read data
if exist(pdfile, 'file')
    fprintf('Read rawdata from: %s\n', pdfile);
    [dataflow, protocol] = CRIS2dataflow(pdfile, CRISpath);
else
    fprintf(2, 'Can not find the rawdata file: %s\n', pdfile);
    return;
end

% output file name
rawdatafile = fullfile(outpath, [pdname, '_' versionflag '.raw']);
fprintf('Transfer the rawdata to: %s\n', rawdatafile);

% fix offset
Z0 = 16384;
Raw_Data = dataflow.rawdata - mean(dataflow.offset, 2) + Z0;
Raw_Data = num2cell(Raw_Data, 1);

% values to put in rawdata 
versionstr = regexp(versionflag, '(\.\d+)|(v\d+)', 'match');
rawdataversion = [str2double(versionstr{1}(2:end)), str2double(versionstr{2}(2:end))];
statusflag = hex2dec('8000');
Reading_Number = num2cell(dataflow.rawhead.Reading_Number, 1);
Angle_encoder = num2cell(dataflow.rawhead.Angle_encoder, 1);
Integration_Time = num2cell(dataflow.rawhead.Integration_Time, 1);
KV = num2cell(dataflow.rawhead.KV, 1);
mA = num2cell(dataflow.rawhead.mA, 1);
Nview = protocol.numView;
rawdatasize = size(dataflow.rawdata, 1)*3;    % 24bit

% raw
raw = struct();
raw(Nview) = struct();
[raw(:).Package_Version] = deal(rawdataversion);
[raw(:).Status_Flag] = deal(statusflag);
[raw(:).Series_Number] = deal(1);
[raw(:).Shot_Number] = deal(1);
[raw(:).Reading_Number] = Reading_Number{:};
[raw(:).Angle_encoder] = Angle_encoder{:};
[raw(:).Integration_Time] = Integration_Time{:};
[raw(:).KV] = KV{:};
[raw(:).mA] = mA{:};
% [raw(:).Start_Slice] = deal(startslice); TBC
% [raw(:).End_Slice] = deal(endslice); TBC
% [raw(:).Slice_merge] = deal(slicemerge); TBC
[raw(:).Slice_Number] = deal(protocol.numRow);
[raw(:).Raw_Data_Size] = deal(rawdatasize);
[raw(:).Raw_Data] = Raw_Data{:};

% find the format configure file
rawcfgfile = cfgmatchrule(rawdatafile, cfgpath, versionflag);
rawcfg = readcfgfile(rawcfgfile);
% pack the data (to file)
packstruct(raw, rawcfg, rawdatafile);

end