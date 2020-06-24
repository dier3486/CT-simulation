function rawfile = checkrawfile(rawpath, rawext, allIwant)
% looking for files in the folder rawpath match with ext rawext
% rawfile = checkrawfile(rawpath, rawext, allIwant);
% e.g. 
% rawfile = checkrawfile(rawpath, '.pd'); to return the first .pd file in rawpath
% rawfile = checkrawfile('./data/**', '.(raw|pd)', true); to return all the .raw and .pd files in extension zero or more 
% folders under ./data/

if nargin<2
    rawext = '.raw';
end
if nargin<3
    allIwant = false;
end

if exist(rawpath, 'file') == 2
    % it is a file
    rawfile = rawpath;
    return;
else
    rawfile = getextfilesinfolader(rawpath, rawext);
    if ~allIwant
        if ~isempty(rawfile)
            rawfile = rawfile{1};
        else
            rawfile = '';
        end
    end
end

end


function filename = getextfilesinfolader(filepath, fileext)

pathdir = dir(filepath);
fileextmatch = ['.*\' fileext '$'];
filename = regexpi({pathdir.name}, fileextmatch, 'match');
% filename = [filename{:}];
filename = fullfile({pathdir(~cellfun(@isempty, filename)).folder}, [filename{:}]);
% NOTE: when the fileext is '.raw' the fileextmatch is '.*\.raw$' which is equavalent with dir *.raw;
% and to use the fileext='.(raw|pd)' to do dir *.raw and *.pd. 
end
