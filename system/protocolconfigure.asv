function protocol = protocolconfigure(protocol_cfg)
% configure protocol

% merge
protocol = structmerge(emptyprotocol(), protocol_cfg, 1);

% scan
if isempty(protocol.scan)
    % do not to scan
    1;
    % do nothing
end
% collimator
if isempty(protocol.collimator)
    % collimator is close
    1;
    % do nothing
end
% bowtie
if isempty(protocol.bowtie)
    % empty bowtie
    protocol.bowtie = 0;
end
% focalspot
if isempty(protocol.focalspot)
    % default focalspot
    protocol.focalspot = 1;
end
% focalsize
if isempty(protocol.focalsize)
    % default focalsize
    protocol.focalsize = 1;
end
% KV mA
if isempty(protocol.KV)
    % no Xray
    protocol.KV = 0;
end
if isempty(protocol.mA)
    % no Xray
    protocol.mA = 0;
end
% multi-KVmA
protocol.KV = protocol.KV(:)';
protocol.mA = protocol.mA(:)';
% double check size
lengthKV = length(protocol.KV);
lengthmA = length(protocol.mA);
if lengthKV>1 && lengthmA>1 && lengthKV~=lengthmA
    error('KV mA number is not match in protocol.');
end
% mA_air
if isempty(protocol.mA_air)
    % default air scan
    protocol.mA_air = protocol.mA;
end


end


function protocol = emptyprotocol()
% what a protocol shall include
protocol = struct();

protocol.scan = [];
protocol.collimator = [];
protocol.bowtie = [];
protocol.focalspot = [];
protocol.focalsize = [];
protocol.KV = [];
protocol.mA = [];
protocol.mA_air = [];
protocol.viewperrot = [];
protocol.rotationspeed = [];
protocol.rotationnumber = [];
protocol.viewnumber = [];
protocol.startangle = [];
protocol.startcouch = [];
protocol.shotnumber = [];
protocol.shotcouchstep = [];
protocol.couchheight = [];
protocol.couchspeed = [];
protocol.rawdatastyle = [];

end