function R = xml2json(filename, recurse_flag, overwrite_flag)

if nargin < 3
    overwrite_flag = true;
end
if nargin < 2
    recurse_flag = true;
end

if isstring(filename)
    filename = char(filename);
end

R = 0;
if isfile(filename)
    [filepath,name,ext] = fileparts(filename);
    if strcmpi(ext, '.xml')
        jsonfile = fullfile(filepath, [name, '.json']);
        if ~isfile(jsonfile) || overwrite_flag
            fprintf('%s  -->  .json\n', filename);
            S = myxml2struct(filename);
            jsonwrite(S, jsonfile);
            R = R + 1;
        end
    end
elseif isfolder(filename)
    filedir = dir(filename);
    for ii = 1:length(filedir)
        filename_ii = fullfile(filedir(ii).folder, filedir(ii).name);
        if filedir(ii).isdir
            if ~any(strcmp(filedir(ii).name, {'.', '..'})) && recurse_flag
                R1 = xml2json(filename_ii, recurse_flag, overwrite_flag);
            else
                R1 = 0;
            end
        else
            R1 = xml2json(filename_ii, recurse_flag, overwrite_flag);
        end
        R = R + R1;
    end
end

end