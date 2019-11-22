function [samplekeV, viewangle, couch] = scanprepare(SYS)
% scan prepare, called in projection scan

% parameters to use
Nview_pr = SYS.protocol.viewperrot;
Nview = SYS.protocol.viewnumber;
Nrot = SYS.protocol.rotationnumber;
Nshot = SYS.protocol.shotnumber;
startangle = SYS.protocol.startangle;

% samplekeV
if strcmpi(SYS.simulation.spectrum, 'Single')
    samplekeV = SYS.world.refrencekeV;
else
    samplekeV = SYS.world.samplekeV;
end

% viewangles
switch lower(SYS.protocol.scan)
    case {'axial', 'helical'}
        % rotation
        viewangle = 0 : pi*2/Nview_pr : pi*2*Nrot;
        viewangle = viewangle(1:Nview);
    otherwise
        % no rotation
        viewangle = zeros(1, Nview);
end
% multi shot
viewangle = reshape(repmat(viewangle(:), 1, Nshot) + startangle, 1, []);
% viewangle = mod(viewangle, pi*2);

% couch (TBC)
couch = 0;

end