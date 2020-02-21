% load raw 

% addpath(genpath('../CRIS'));
curpwd = pwd;

% rawpath = 'E:\data\rawdata\spect_data\';
rawpath = 'F:\data-Dier.Z\PG\Cu\';
dirpath = dir(rawpath);

% Data files
Data = struct();
index = 1;
for ii=1:length(dirpath(:))
    path_ii = dirpath(ii);
    if ~path_ii.isdir
        continue
    end
    if strcmp(path_ii.name, '.') || strcmp(path_ii.name, '..')
        continue
    end
    subdir = dir(fullfile(path_ii.folder, path_ii.name, '*.pd'));
    if isempty(subdir)
        continue
    end
    Data(index).filename = fullfile(subdir(1).folder, subdir(1).name);
    namekey = regexp(path_ii.name, '_', 'split');
    % I know the keys are:
    Data(index).object = namekey{1};
    KV = regexp(namekey{2}, '\d+(kV|KV)','match');
    Data(index).KV = str2double(KV{1}(1:end-2));
    mA = regexp(namekey{2}, '\d+mA','match');
    Data(index).mA = str2double(mA{1}(1:end-2));
%     Data(index).bowtie = namekey{3};
    Data(index).focal = namekey{3};
    % 2019/11/21
    index = index + 1;
end

cd('../CRIS');
% load data
for ii = 1:length(Data)
    fprintf('Reading %s ...', Data(ii).filename);
    raw_ii = readRaw(Data(ii).filename, 1);
    Data(ii).offsetmean = squeeze(mean(raw_ii.offset.chanData, 1));
    Nskip = 200;
    Data(ii).rawmean = squeeze(mean(raw_ii.active.chanData(Nskip:end,:,:), 1));
    Data(ii).viewmean = squeeze(mean(raw_ii.active.chanData, 2));
    fprintf('done\n');
end
cd(curpwd);