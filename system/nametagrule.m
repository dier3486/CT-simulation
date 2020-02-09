function nametags = nametagrule(namerule, protocol, KV, mA)
% file name rule, usually for .raw, .corr and reconxml
% nametags = nametagrule(namerule, protocol, KV, mA);
% or nametags = nametagrule(namerule, protocol);
% Then a typical rawdata filename is like this: ['rawdata' namekey nametags '_' 'v1.0' '.raw'],
% INPUT
%   namerule        to select a name rule, 'standard', 'simple', 'timestamp' or otherwise as default
%   protocol        protocol struct of simulation or recon which could be SYS.protocol or reconxml.protocol
%   KV, mA          in simulation we might once run multi KVmAs, in that case we can loop the function to get multi file names 

if nargin > 3
    KVmA = ['_' num2str(KV) 'KV' num2str(mA) 'mA'];
elseif nargin > 2
    KVmA = ['_' num2str(KV) 'KV'];
else
    KVmA = '';
end

switch lower(namerule)
    case 'standard'
        % standard name rule
        nametags = ['_' protocol.scan '_' protocol.bowtie '_' protocol.collimator ...
                KVmA '_' num2str(protocol.rotationspeed) 'secprot'];
    case 'simple'
        % only series number
        if length(protocol.KV)>1 || length(protocol.mA)>1
            nametags = ['_series' num2str(protocol.series_index) KVmA];
        else
            nametags = ['_series' num2str(protocol.series_index)];
        end
    case {'time', 'timestamp'}
        % time stamp
        nametags = num2str(now, '%.10f');
        pause(0.001);
    otherwise
        % series number and KVmA
        nametags = ['_series' num2str(protocol.series_index) KVmA];
end

end