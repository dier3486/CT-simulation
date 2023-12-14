function [samplekeV, viewangle, couch, shotindex, gantrytilt] = scanprepare(SYS, iblock)
% scan prepare, called in projection scan

% parameters to use
Nview_pr = SYS.protocol.viewperrot;
Nview = SYS.protocol.viewnumber;    % Nview per shot
Nrot = SYS.protocol.rotationnumber;
Nshot = SYS.protocol.shotnumber;
startangle = SYS.protocol.startangle;
rotdirect = SYS.protocol.rotationdirect;
startcouch = SYS.protocol.startcouch;
couchstep = SYS.protocol.shotcouchstep;
couchspeed = SYS.protocol.couchspeed;
couchheight = SYS.protocol.couchheight;
rotspeed = SYS.protocol.rotationspeed;
inttime = SYS.protocol.integrationtime;
gantrytilt = SYS.protocol.gantrytilt;
Startview = SYS.console.Startviewperblk(iblock);
if iblock < SYS.console.Nviewblk
    Nviewofcurrblk = SYS.console.viewblock;
else
    Nviewofcurrblk = Nview * Nshot - SYS.console.viewblock * (SYS.console.Nviewblk-1);
end

% samplekeV
if strcmpi(SYS.simulation.spectrum, 'Single')
    samplekeV = SYS.world.referencekeV;
else
    samplekeV = SYS.world.samplekeV;
end

% startangle to pi
startangle = mod(startangle*(pi/180), pi*2);

% viewangles
switch lower(SYS.protocol.scan)
    case {'axial', 'helical'}
        % rotation
        viewangle = (0 : pi*2/Nview_pr : pi*2*Nrot).*rotdirect;
        viewangle = viewangle(1:Nview);
    otherwise
        % no rotation
        viewangle = zeros(1, Nview);
end
% multi shot + startangle
viewangle = reshape(repmat(viewangle(:), 1, Nshot) + startangle, 1, []);
% viewangle = mod(viewangle, pi*2);
% block
viewangle = viewangle(Startview : Startview + Nviewofcurrblk-1);


% gantry tilt
% to pi
gantrytilt = gantrytilt(:).*(pi/180);
% multi shot
if size(gantrytilt, 1) == 1
    gantrytilt = repmat(gantrytilt, Nshot, 1);
end
% multi shot of tilt
gantrytilt = repelem(gantrytilt, Nviewofcurrblk, 1);

% shotindex
shotindex = reshape(repmat(1:Nshot, Nview, 1), 1, []);
% block
shotindex = shotindex(Startview : Startview + Nviewofcurrblk-1);

% couch
switch lower(SYS.protocol.scan)
    case 'axial'
        couch_z = repmat((0:Nshot-1).*couchstep + startcouch, Nview, 1);
    case 'helical'
        couch_z = (1:Nview)'.*(couchspeed*rotspeed/Nview_pr) + ...
            (0:Nshot-1).*couchstep + startcouch;
    otherwise
        % static or topo
        couch_z = (1:Nview)'.*(couchspeed*inttime*1e-6) + ...
            (0:Nshot-1).*couchstep + startcouch;
end
couch = [zeros(Nview*Nshot, 1) -ones(Nview*Nshot, 1).*couchheight couch_z(:)];
% block
couch = couch(Startview : Startview + Nviewofcurrblk-1, :);

end