function corrfile = corrcouplerule(protocol, corrpath, corrcouple, corrname, fileext)
% looking for corr table
% corrfile = corrcouplerule(protocol, corrpath, corrcouple, corrname, fileext);
% to find the file in folder corrpath which matched with the protocol by the match rule corrcouple.(corrname).
% INPUT
%   protocol                the struct of protocol, e.g. the reconxml.recon.protocol which configured buy a recon .xml.
%    

if nargin<5
    % default fileext
    fileext = '.corr';
end

% load corr couple rule configure file
if ~isstruct(corrcouple)
    corrcouple = readcfgfile(corrcouple);
end
if isfield(corrcouple, corrname)
    corrtags = fieldnames(corrcouple.(corrname));
    Ntag = length(corrtags);
else
    error('Not defined corr couple rule for %s!', corrname);
end
% get tagnames to couple
tagnames = corrtags;
% default tagnames is the corrtags
for itag = 1:Ntag
    % but user can set it to another name
    if isfield(corrcouple.(corrname).(corrtags{itag}), 'tagname')
        tagnames{itag} = corrcouple.(corrname).(corrtags{itag}).tagname;
    end
end
% get the tag strings from protocol of the tagnames
[protocoltags, matchfun] = protocol2nametag(protocol, tagnames);
% we will check if the strings in protocoltags{} exist in the filenames by the matchfun{}.
% A typical matchfun is @strcmp, and user can set it to any function_handle.

% filename is dir *.fileext
pathdir = dir(corrpath);
fileextmatch = ['.*\' fileext '$'];
filename = regexpi({pathdir.name}, fileextmatch, 'match');
filename = [filename{:}];
% NOTE: when the fileext is '.corr' the fileextmatch is '.*\.corr$' which is to do dir *.corr;
% and to use the fileext='.(corr|ct)' to do dir *.corr and *.ct. 

% loop the filename to find out which is (most) match with the protocol
Nfile = size(filename(:), 1);
filecouple = zeros(Nfile, 1);
for ifile = 1:Nfile
    % check if the ifile-th file matched
    % split the filename
    [~, filename_part, ~] = fileparts(filename{ifile});
    filenametags = regexp(filename_part, '_', 'split');
    % corrname
    if ~strcmpi(filenametags{1}, corrname) && ~isfield(corrcouple.(corrname), 'corrname')
        filecouple(ifile) = -1;
        continue;
    end
    for itag = 1:Ntag
        % check the couple value of the protocoltags{itag}
        couplevar = tagcouple(filenametags, protocoltags{itag}, corrcouple.(corrname).(corrtags{itag}), matchfun{itag});
        % -1: not couple; 0: perfect couple; >0 toleranced couple, less is better
        if couplevar<0
            % not couple
            filecouple(ifile) = -1;
            break;
        else
            filecouple(ifile) = filecouple(ifile) + couplevar;
        end
    end
    % if all the protocoltags are coupled, the file is matched
end

% is any file matched?
s = find(filecouple>=0);
if ~isempty(s)
    % return the (most) matched file
    [~, tmin] = min(filecouple(s));
    corrfile = fullfile(corrpath, filename{s(tmin)});
else
    corrfile = '';
end

end


function couplevar = tagcouple(filenametags, protocoltag, ruletag, matchfun)

% split to sopport multi protocoltag
protocoltag_split = regexp(protocoltag, '_', 'split');
Nsplit = length(protocoltag_split);
if Nsplit == 1
    % call tagcouple2 to judge if the protocoltag match with the filenametags
    couplevar = tagcouple2(filenametags, protocoltag, ruletag, matchfun);
else  % Nsplit>1
    % to support the protocoltag with '_', e.g. 'water300_center'
    couplevar = zeros(1, Nsplit);
    for itag = 1:Nsplit
        couplevar(itag) = tagcouple2(filenametags, protocoltag_split{itag}, ruletag, matchfun);
    end
    if any(couplevar==-1)
        couplevar = -1;
    else
        couplevar = sum(couplevar);
    end
end

end


function couplevar = tagcouple2(filenametags, protocoltag, ruletag, matchfun)
% subfunction to check if the nametags match with protocoltag 

% default matchfun
if nargin<4
    matchfun = @strcmpi;
end
% replace the matchfun if it was defined in ruletag
if isfield(ruletag, 'matchfun') && isa(ruletag.matchfun, 'function_handle')
    matchfun = ruletag.matchfun;
end
% NOTE: the matchfun could be @strcmp for a full compare, @strcmpi to ignore case or a function_handle, e.g. to call regexp to
% define more flexible comparing rules, 
% e.g. to set matchfun = @(x, y) ~isempty(regexpi(x, y)) as a partly matched ignoring case rule.

% check if the nametags match with protocoltag
if any(cellfun(@(x) matchfun(x, protocoltag), filenametags))
    % matched perfectly
    couplevar = 0;
elseif isfield(ruletag, 'match')
    % seems not exactly matched, but try a tolerenced match
    if ~iscell(ruletag.match)
        toledmatchrule = {ruletag.match};
    else
        toledmatchrule = ruletag.match;
    end
    matchwith = cellfun(@(x) tagintolmatch(filenametags, protocoltag, x, matchfun), toledmatchrule);
    s = matchwith>0;
    if any(s)
        % matched in tolerenced level min(matchwith(s))
        couplevar = min(matchwith(s));
    else
        % not match even with tolerence
        couplevar = -1;
    end
else
    % the tag is not match
    couplevar = -1;
end

end


function r = tagintolmatch(filenametags, protocoltag, ruletag, matchfun)
% subfunction to check if the nametags match with protocoltag in tolerence ruletag.matchwith
if isempty(ruletag)
    r = -1;
    return;
end

% split matchwith tags
if ~iscell(ruletag.matchwith)
	matchwithtags = regexp(ruletag.matchwith, '(, +)|(,)', 'split');
else
    matchwithtags = ruletag.matchwith;
end

if isfield(ruletag, 'whichin')
    % split whichin tags
    if ~iscell(ruletag.whichin)
        whichintags = regexp(ruletag.whichin, '(, +)|(,)', 'split');
    else
        whichintags = ruletag.whichin;
    end
    % if the protocoltag in whichintags
    if any(cellfun(@(x) strcmpi(x, protocoltag), whichintags))
        % if the matchwithtag in nametags
        r = -1;
        for ii = 1:length(matchwithtags)
            if any(cellfun(@(x) matchfun(x, matchwithtags{ii}), filenametags))
                % matched in tolerence level ii
                r = ii;
                break
            end
        end
    else
        % not match
        r = -1;
    end
else
    % if the matchwithtag in nametags
    r = -1;
    for ii = 1:length(matchwithtags)
        if any(cellfun(@(x) matchfun(x, matchwithtags{ii}), filenametags))
            % matched in tolerence level ii
            r = ii;
            break
        end
    end
end

end

