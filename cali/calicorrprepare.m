function datafile = calicorrprepare(datafile, filepath, fileext, newdatafield)
% prepare .corr for calibration
% datafile = calicorrprepare(datafile, filepath);
% or, datafile = calicorrprepare(datafile, filepath, fileext, newdatafield);
% sorry for no more helps, plz look up the cali scripts for more information.

if nargin<3
    % default fileext
    fileext = '.corr';
end
if nargin<4
    newdatafield = 'calitable';
end


% size of datafile
Ndata = size(datafile(:), 1);

% toloop tags
loopfields = fieldnames(datafile);
Nloop = size(loopfields(:),1);

% files
filefields = fieldnames(filepath);
Nfield = size(filefields(:), 1);

for ifield = 1:Nfield
    dataname = filefields{ifield};
    pathdir = dir(filepath.(dataname).path);
    Npd = size(pathdir(:), 1);
    % namekey
    if isfield(filepath.(dataname), 'namekey') && ~isempty(filepath.(dataname).namekey)
        namekey = filepath.(dataname).namekey;
        if ~iscell(namekey)
            namekey = {namekey};
        end
    else
        namekey = {};
    end
    % skip
    if isfield(filepath.(dataname), 'skiptag')  && ~isempty(filepath.(dataname).skiptag)
        skiptag = filepath.(dataname).skiptag;
        if ~iscell(skiptag)
            skiptag = {skiptag};
        end
    else
        skiptag = {};
    end
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
            if ~isempty(namekey)
                % is all the namekey(s) in nametags?
                isnamekeyintags = all(cellfun(@(y) any(cellfun(@(x) strcmpi(x, y), nametags)), namekey));
                if ~isnamekeyintags
                    % namekey not match
                    continue
                end
            end
            % check the tags in toloop
            ismatch = true;
            for itag = 1:Nloop
                tagtoloop = loopfields{itag};
                % is empty?
                if isempty(datafile(ii).(tagtoloop))
                    % empty tag == match
                    continue;
                end
                % is struct?
                if isstruct(datafile(ii).(tagtoloop))
                    % struct == match
                    continue;
                end
                % is skip?
                if ~isempty(skiptag)
                    % is the tagtoloop in skiptag(s)?
                    isskiptag = any(cellfun(@(x) strcmpi(x, tagtoloop), skiptag));
                    if isskiptag
                        % skip == match
                        continue
                    end
                end
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
                datafile(ii).(newdatafield).(dataname) = filename;
            end
        end
    end
end

end