function datafile = calidataprepare(toloop, filepath, fileext)
% prepare data for calibration

if nargin<3
    % fileext
    fileext = '.raw';
end

% expand the values in toloop to datafile
datafile = loopcurse(toloop);
Ndata = size(datafile(:), 1);

% toloop tags
loopfields = fieldnames(toloop);
Nloop = size(loopfields(:),1);

% files
filefields = fieldnames(filepath);
Nfield = size(filefields(:), 1);

for ifield = 1:Nfield
    dataname = filefields{ifield};
    pathdir = dir(filepath.(dataname).path);
    Npd = size(pathdir(:), 1);
    for ipd = 1:Npd
        dirname = pathdir(ipd).name;
        % skip . and ..
        if strcmp(dirname, '.') || strcmp(dirname, '..')
            continue;
        end
        [~, dirpartname, dirext] = fileparts(dirname);
        if pathdir(ipd).isdir
            % dir again
            diragain = dir(fullfile(pathdir(ipd).folder, dirname, ['*' fileext]));
            if isempty(diragain)
                continue;
            else
                filename = fullfile(diragain(1).folder, diragain(1).name);
            end
        elseif strcmpi(dirext, fileext)
            filename = fullfile(pathdir(ipd).folder, dirname);
            dirname = dirpartname;
        else
            continue;
        end
        nametags = regexp(dirname, '_', 'split');
        for ii = 1:Ndata
            % check namekey
            if isfield(filepath.(dataname), 'namekey') && ~isempty(filepath.(dataname).namekey)
                namekey = filepath.(dataname).namekey;
                if ~any(cellfun(@(x) strcmpi(x, namekey), nametags))
                    % namekey not match
                    continue;
                end
            end
            % check the tags in toloop
            ismatch = true;
            for itag = 1:Nloop
                tagtoloop = loopfields{itag};
                switch tagtoloop
                    case 'KV'
                        % check KV
                        KVtag = [num2str(datafile(ii).KV) 'KV'];
                        KVmatch = ~cellfun(@isempty, regexp(nametags, KVtag));
                        if ~any(KVmatch)
                            % KV not match
                            ismatch = false;
                            continue;
                        end
                    case 'focalsize'
                        % check focal size
                        focaltag = [datafile(ii).focalsize 'Focal'];
                        focalmatch = ~cellfun(@isempty, regexpi(nametags, focaltag));
                        if ~any(focalmatch)
                            % focalsize not match
                            ismatch = false;
                            continue;
                        end
                    otherwise
                        % check other tags (bowtie, collimator, ...)
                        checktag = datafile(ii).(tagtoloop);
                        if ~any(cellfun(@(x) strcmpi(x, checktag), nametags))
                            % the tag is not match
                            ismatch = false;
                            continue;
                        end
                end
            end
            if ismatch
                % data matched
                datafile(ii).filename.(dataname) = filename;
            end
        end
    end
end

end


function [datafile, idata] = loopcurse(toloop, datafile, ifield, idata)
% expand the values in toloop to datafile

if nargin<2
    datafile = struct();
end
if nargin<3
    ifield = 1;
    idata = 1;
end

toloopfields = fieldnames(toloop);
Nfields = size(toloopfields(:), 1);

currfield = toloop.(toloopfields{ifield});
Nlist = size(currfield(:), 1);
if Nlist==1 && ~iscell(currfield)
    currfield = {currfield};
end
% loop the elements in current field
for ii = 1:Nlist
    % copy the values in current field to datafile
    if iscell(currfield)
        datafile(idata).(toloopfields{ifield}) = currfield{ii};
    else
        datafile(idata).(toloopfields{ifield}) = currfield(ii);
    end
    % to recurse
    if ifield<Nfields
        [datafile, idata] = loopcurse(toloop, datafile, ifield+1, idata);
    end
    % next
    if ii<Nlist
        datafile(idata+1) = datafile(idata);
        idata = idata+1;
    end
end

end