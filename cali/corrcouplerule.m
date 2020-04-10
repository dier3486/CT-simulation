function corrfile = corrcouplerule(protocol, corrpath, corrcouple, corrname, fileext)
% looking for corr table
% corrfile = corrcouplerule(protocol, corrpath, corrcouple, corrname, fileext);
% to find to the file name ('corrfile') of a corrtable of 'corrname' by 'protocol' in path 'corrpath' with the couple rule
% 'corrcouple'.
% INPUT
%   protocol                the struct of protocol, e.g. the reconxml.recon.protocol which configured buy a recon .xml.
%    

if nargin<5
    % default fileext
    fileext = '.corr';
end

% load corr couple rule configure file
if ischar(corrcouple)
    corrcouple = readcfgfile(corrcouple);
end
if isfield(corrcouple, corrname)
    tagnames = fieldnames(corrcouple.(corrname));
    Ntag = size(tagnames(:), 1);
else
    error('Not defined corr couple rule for %s!', corrname);
end

pathdir = dir(fullfile(corrpath, ['*' fileext]));
Nfile = size(pathdir(:), 1);

filecouple = zeros(Nfile, 1);
for ifile = 1:Nfile
    % check if the file coupled
    dirname = pathdir(ifile).name;
    nametags = regexp(dirname, '_', 'split');
    % corrname
    if ~strcmpi(nametags{1}, corrname)
        filecouple(ifile) = -1;
        continue;
    end
    for itag = 1:Ntag
        if ~isfield(protocol, tagnames{itag})
            continue;
        end
        
        switch tagnames{itag}
            case 'focalspot'
                switch protocol.focalspot
                    case 1
                        protocoltag = 'FocalQFS';
                    case 2
                        protocoltag = 'FocalQFS';
                    case 3
                        protocoltag = 'FocalDFS';
                    otherwise
                        protocoltag = ['Focal' protocol.focalspot];
                end
                matchfun = @(x, y) ~isempty(regexpi(x, y));
            case 'focalsize'
                switch protocol.focalsize
                    case 1
                        protocoltag = 'smallFocal';
                    case 2
                        protocoltag = 'largeFocal';
                    otherwise
                        protocoltag = [protocol.focalsize 'Focal'];
                end
                matchfun = @(x, y) ~isempty(regexpi(x, y));
            case 'KV'
                protocoltag = [num2str(protocol.KV) 'KV'];
                matchfun = @(x, y) ~isempty(regexp(x, y));
            case 'mA'
                protocoltag = [num2str(protocol.mA) 'mA'];
                matchfun = @(x, y) ~isempty(regexp(x, y));
            case 'rotationspeed'
                protocoltag = [num2str(protocol.rotationspeed) 'SecpRot'];
                matchfun = @strcmpi;
            otherwise
                protocoltag = protocol.(tagnames{itag});
                matchfun = @strcmpi;
        end
        couplevar = tagcouple(nametags, protocoltag, corrcouple.(corrname).(tagnames{itag}), matchfun);
        if couplevar<0
            % not match
            filecouple(ifile) = -1;
            break;
        else
            filecouple(ifile) = filecouple(ifile) + couplevar;
        end
    end
    
end

s = find(filecouple>=0);
if ~isempty(s)
    [~, tmin] = min(filecouple(s));
    corrfile = fullfile(corrpath, pathdir(s(tmin)).name);
else
    corrfile = '';
end

end


function couplevar = tagcouple(nametags, protocoltag, ruletag, matchfun)
% subfunction to check if the nametags match with protocoltag 
if nargin<4
    matchfun = @strcmpi;
end

% check if the nametags match with protocoltag
if any(cellfun(@(x) matchfun(x, protocoltag), nametags))
    % matched perfectly
    couplevar = 0;
elseif ~isempty(ruletag)
    % seems not exactly matched, but try a tolerenced match
    if ~iscell(ruletag)
        ruletag = {ruletag};
    end
    matchwith = cellfun(@(x) tagintolmatch(nametags, protocoltag, x, matchfun), ruletag);
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


function r = tagintolmatch(nametags, protocoltag, ruletag, matchfun)
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
            if any(cellfun(@(x) matchfun(x, matchwithtags{ii}), nametags))
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
        if any(cellfun(@(x) matchfun(x, matchwithtags{ii}), nametags))
            % matched in tolerence level ii
            r = ii;
            break
        end
    end
end

end

