datapath = '\\192.168.1.36\pangu.dat\RAW\smallbowtie';

fdir = dir(fullfile(datapath, '**\*.pd'));
Nf = size(fdir(:) ,1);
% Nf = 2;

% angle sections
Nsection = 12;
delta = pi*2/Nsection;
sectangle = (0:Nsection).*delta;

rawdata = cell(Nf, 1);
for ii = 1:Nf
    fprintf('.');
    dataflow = CRIS2dataflow(fullfile(fdir(ii).folder, fdir(ii).name));
    % viewangle
    viewangle = dataflow.rawhead.viewangle;
    viewangle = mod(viewangle + delta/2, pi*2);
    % ini rawdata
    rawdata{ii} = zeros(size(dataflow.rawdata, 1), Nsection);
    for isect = 1:Nsection 
        viewindex = (viewangle>=sectangle(isect)) & (viewangle<sectangle(isect+1));
        rawdata{ii}(:, isect) = reshape(mean(dataflow.rawdata(:, viewindex), 2), [], 1);
    end
end

fprintf('\n');