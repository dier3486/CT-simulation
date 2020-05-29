function rawfile = checkrawfile(rawfile, rawext, allIwant)
% a sub function to look for rawfile in a folder
% rawfile = checkrawfile(rawfile, rawext, allIwant);
% or, rawfile = checkrawfile(rawfile, '.pd');

if nargin<2
    rawext = '.raw';
end
if nargin<3
    allIwant = false;
end

if exist(rawfile, 'file') == 2
    % do nothing
    return;
elseif exist(rawfile, 'dir')
    filename = getextfilesinfolader(rawfile, rawext);
    if ~isempty(filename)
        if allIwant
            rawfile = fullfile(rawfile, filename);
        else
            rawfile = fullfile(rawfile, filename{1});
        end
    else
        rawfile = '';
    end
else
    rawfile = '';
end

end


function filename = getextfilesinfolader(filepath, fileext)

pathdir = dir(filepath);
fileextmatch = ['.*\' fileext '$'];
filename = regexpi({pathdir.name}, fileextmatch, 'match');
filename = [filename{:}];
% NOTE: when the fileext is '.raw' the fileextmatch is '.*\.raw$' which is equavalent with dir *.raw;
% and to use the fileext='.(raw|pd)' to do dir *.raw and *.pd. 
end
