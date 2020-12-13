% load raw 

rawpath = 'D:\data\taiying\CU\';
dirpath = dir(fullfile(rawpath, '*.raw'));

% Data files
Data = struct();
index = 1;
for ii=1:length(dirpath(:))
    
    Data(index).filename = fullfile(dirpath(ii).folder, dirpath(ii).name);
    [~, filename_tags, ~] = fileparts(Data(index).filename);
    namekey = regexp(filename_tags, '_', 'split');
    % I know the keys are:
    Data(index).object = namekey{end};
    KV = regexp(namekey, '\d+(kv|KV)','match');
    KV = KV{~cellfun(@isempty, KV)};
    Data(index).KV = str2double(KV{1}(1:end-2));
    mA = regexp(namekey, '\d+ma','match');
    mA = mA{~cellfun(@isempty, mA)};
    Data(index).mA = str2double(mA{1}(1:end-2));
%     Data(index).bowtie = namekey{3};
%     Data(index).focal = namekey{3};
    % 2019/11/21
    index = index + 1;
end

% load data
rawcfgfile = 'D:\matlab\ct\BCT16\IOstandard\rawdata_v1.0.xml';
for ii = 1:length(Data)
    fprintf('Reading %s ...', Data(ii).filename);
    raw_ii = loadbindata(Data(ii).filename, rawcfgfile);
    Nskip = 100;
    Data(ii).rawmean = mean([raw_ii(Nskip+1:end).Raw_Data], 2);
    fprintf('done\n');
end
