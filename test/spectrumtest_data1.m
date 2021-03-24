path1 = 'F:\data-Dier.Z\NV\ÄÜÆ×ÏìÓ¦²âÊÔ\dark(offset)';

dir1 = dir(path1);
% Nthick = size(dir1,1)-2;
% thickness = zeros(Nthick, 1);
thickness = [];
data = [];
for ii = 3:size(dir1,1)
    if dir1(ii).isdir
        % take thickness
        tokenmm = regexp(dir1(ii).name, '(.+)mm', 'tokens');
        if ~isempty(tokenmm)
            thickness = [thickness str2double(tokenmm{1}{1})];
        end
        % load data
        datafiles = dir(fullfile(path1, [dir1(ii).name '/*.raw']));
        data_ii = 0;
        for jj = 1:size(datafiles, 1)
            fid = fopen(fullfile(datafiles(jj).folder, datafiles(jj).name), 'r');
            data_jj = fread(fid, inf, 'uint16');
            fclose(fid);
            data_ii = data_ii + double(data_jj);
        end
        data_ii = data_ii./size(datafiles, 1);
        data = [data data_ii(:)];
    end
end