function [prototag, matchfun] = protocol2nametag(protocol, tags)
% get tags from protocol
% [prototag, matchfun] = protocol2nametag(protocol, tags);
% or
% prototag = protocol2nametag(protocol, tags);

if nargin<2
    tags = [];
end

% % default tags (NO)
% default_tags = {'series', 'scan', 'bowtie', 'collimator', 'KV', 'mA', 'Focalsize', 'Focalspot', 'rotationspeed'};
% if isempty(tags)
%     tags = default_tags;
% end

Ntag = size(tags(:),1);
prototag = cell(1, Ntag);
matchfun = cell(1, Ntag);

for itag = 1:Ntag
    % default matchfun
    matchfun{itag} = @strcmp;
    % switch the tags
    switch lower(tags{itag})
        case {'series', 'seriesindex'}
            prototag{itag} = ['series' num2str(protocol.seriesindex)];
        case 'focalsize'
            prototag{itag} = tagfocalsize(protocol.focalsize, tags{itag});
            matchfun{itag} = @(x, y) ~isempty(regexpi(x, y));
        case 'focalspot'
            prototag{itag} = tagfocalspot(protocol.focalspot);
            matchfun{itag} = @(x, y) ~isempty(regexp(x, y, 'once'));
        case {'excbowtie', 'largebowtie'}
            % escape bowtie to air, large and small
            prototag{itag} = exclargebowtie(protocol.bowtie);
        case 'kv'
            prototag{itag} = [num2str(protocol.KV) tags{itag}];
            matchfun{itag} = @(x, y) ~isempty(regexp(x, y, 'once'));
        case 'ma'
            prototag{itag} = [num2str(protocol.mA) tags{itag}];
            matchfun{itag} = @(x, y) ~isempty(regexp(x, y, 'once'));
        case 'kvma'
            % put KV mA together
            prototag{itag} = [num2str(protocol.KV) 'KV' num2str(protocol.mA) 'mA'];
        case 'focal'
            % put focalsize focalspot together
            prototag{itag} = [tagfocalsize(protocol.focalsize) tagfocalspot(protocol.focalspot)];
        case 'rotsec'
            prototag{itag} = [sprintf('%0.2f',protocol.rotationspeed) 'sec'];
        case {'time', 'timestamp'}
            prototag{itag} = num2str(now, '%.10f');
            pause(0.001);
        case 'timeseries'
            % series number . time stamp 
            prototag{itag} = [num2str(protocol.seriesindex) '.' num2str(now, '%.10f')];
            pause(0.001);
        otherwise
            % try to define the tag by general rule
            if isfield(protocol, tags{itag})
                % scan, bowtie, collimator or other
                if ischar(protocol.(tags{itag}))
                    prototag{itag} = protocol.(tags{itag});
                else
                    prototag{itag} = num2str(protocol.(tags{itag}));
                end
            else
                tagsplit = regexp(tags{itag}, '_', 'split');
                if isfield(protocol, tagsplit{1})
                    prototag{itag} = sprintf(tagsplit{2}, protocol.(tagsplit{1}));
                else
                    % copy the string
                    prototag{itag} = tags{itag};
                end
            end
    end
    
end

end


function nametag = tagfocalsize(focalsize, tagkey)
% tag of focal size
switch focalsize
    case 1
        nametag = 'Small';
    case 2
        nametag = 'Large';
    otherwise
        nametag = num2str(focalsize);
end

if nargin>1
    % upper
    if all(isstrprop(tagkey, 'upper'))
        % ALL UPPER
        nametag = [upper(nametag) 'FOCAL'];
    elseif all(isstrprop(tagkey, 'lower'))
        % all lower
        nametag = [lower(nametag) 'focal'];
    elseif isstrprop(tagkey(1), 'upper')
        % First Upper
        nametag(1) = upper(nametag(1));
        nametag = [nametag 'Focal'];
    else
        % unknown request
        nametag = [nametag 'Focal'];
    end
else
    nametag = [nametag 'Focal'];
end

end


function focalspot = tagfocalspot(focalspot)
% tag of focal spot
Nspot = size(focalspot(:), 1);
if Nspot>1
    switch Nspot
        case 2
            focalspot = 'DFS';
        otherwise
            focalspot = num2str(focalspot);
    end
else
    switch focalspot
        case {1, 2}
            focalspot = 'QFS';
        case {3, 6}
            focalspot = 'DFS';
        otherwise
            focalspot = num2str(focalspot);
    end
end

end


function nametag = exclargebowtie(bowtie)
% tag of air, large or small bowtie
switch lower(bowtie)
    case {0, 'empty', 'air'}
        nametag = 'AIRBOWTIE';
    case {1, 'body', 'large'}
        nametag = 'LARGEBOWTIE';
    case {2, 'head', 'small'}
        nametag = 'SMALLBOWTIE';
    otherwise
        % do nothing
        1;
end
end
