function nametags = nametagrule(namerule, protocol, tags, KV, mA)
% file name rule, usually for .raw, .corr and reconxml
% nametags = nametagrule(namerule, protocol, tags, KV, mA);
% or nametags = nametagrule(namerule, protocol, tags);
% Then a typical rawdata filename is like this: ['rawdata' namekey nametags '_' 'v1.0' '.raw'],
% INPUT
%   namerule        to select a name rule, 'standard', 'simple', 'timestamp' or otherwise as default
%   protocol        protocol struct of simulation or recon which could be SYS.protocol or reconxml.protocol
%   tags            tags to appear in name, e.g. {'KV', 'mA', 'collimator'}
%   KV, mA          to define the tag 'KV mA', in simulation we might once run multi KVmAs and loop them to set the output
%                   files' name

% clean tags
if nargin<3
    tags = {};
end
if ~iscell(tags) && ~isempty(tags)
    if ischar(tags)
        tags = regexp(tags, '(, +)|(,)', 'split');
    else
        error('Illegal input tags!');
    end
end

% input KV mA? (only used in simulation task for once calculating multi KV mA)
if nargin > 3
    protocol.KV = KV;
end
if nargin > 4
    protocol.mA = mA';
end

switch lower(namerule)
    case 'standard'
        % standard name rule (quick start for newbie)
        tags = {'scan', 'bowtie', 'collimator', 'KVmA', 'Focal', 'rotationspeed_%gSecpRot'};
    case {'manu', 'manutags'}
        % use the input tags (you should use this to manage your CT files)
        1;
        % do nothing
    case 'simple'
        % series number and KVmA (these are other quick start samples)
        tags = {'series', 'KVmA'};
    case 'series'
        % only series number
        tags = {'series'};
    case {'time', 'timestamp'}
        % time stamp
        tags = {'time'};
    case 'timeseries'
        % series number . time stamp 
        tags = {'timeseries'};
    otherwise
        % empty nametag
        tags = {};
        % WARN: files could be overwrited due to repeated file names.
end

% genrate the tags cell by protocol
prototag = protocol2nametag(protocol, tags);
% link the tags with '_'
nametags = linkstringcell(prototag, '_');

end

