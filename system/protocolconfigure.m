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
switch lower(protocol.focalsize)
    case 'small'
        protocol.focalsize = 1;
    case 'large'
        protocol.focalsize = 2;
    otherwise
        % do nothing
        1;
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
% viewperrot
if isempty(protocol.viewperrot)
    if protocol.rotationnumber>0
        protocol.viewperrot = protocol.viewnumber/protocol.rotationnumber;
    else
        protocol.viewperrot = protocol.viewnumber;
    end
end
% rotationspeed
if isempty(protocol.rotationspeed)
    protocol.rotationspeed = 1;
    % which will be used to calculation the integration time, even for
    % static scan
end
% rotationnumber
if isempty(protocol.rotationnumber)
    % default is one rotation
    protocol.rotationnumber = 1;
end
% viewnumber
if isempty(protocol.viewnumber)
    protocol.viewnumber = floor(protocol.viewperrot*protocol.rotationnumber);
end
% integrationtime
if isempty(protocol.integrationtime)
    protocol.integrationtime = 1.0e6/protocol.rotationspeed/protocol.viewperrot;
    if ~isfinite(protocol.integrationtime) || protocol.integrationtime==0
        error('Wrong integration time! Which due to the rotationspeed/viewperrot is an illegal number.');
    end
end

% startangle
if isempty(protocol.startangle)
    protocol.startangle = 0;    % 0 is tube at -pi/2
end
% startcouch
if isempty(protocol.startcouch)
    protocol.startcouch = 0;
end
% shotnumber
if isempty(protocol.shotnumber)
    protocol.shotnumber = 1;
end
% shotcouchstep
if isempty(protocol.shotcouchstep)
    protocol.shotcouchstep = 0;
end
% couchheight
if isempty(protocol.couchheight)
    protocol.couchheight = 0;
end
% couchspeed
if isempty(protocol.couchspeed)
    protocol.couchspeed = 0;
end
% rawdatastyle
if isempty(protocol.rawdatastyle)
    protocol.rawdatastyle = 'uint24';
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
protocol.integrationtime = [];
protocol.startangle = [];
protocol.startcouch = [];
protocol.shotnumber = [];
protocol.shotcouchstep = [];
protocol.couchheight = [];
protocol.couchspeed = [];
protocol.rawdatastyle = [];

end