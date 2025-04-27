function [samplekeV, viewangle, couch, shotindex, gantrytilt, turningpoints] = scanprepare(SYS, iblock)
% scan prepare, called in projection scan

% parameters to use
Nviewprot = SYS.protocol.viewperrot;
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
    case {'axial', 'helical', 'cradle'}
        % rotation
        viewangle = (0 : pi*2/Nviewprot : pi*2*Nrot).*rotdirect;
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
turningpoints = [];
switch lower(SYS.protocol.scan)
    case 'axial'
        couch_z = repmat((0:Nshot-1).*couchstep + startcouch, Nview, 1);
    case 'helical'
        couch_z = (0:Nview-1)'.*(couchspeed*rotspeed/Nviewprot) + ...
            (0:Nshot-1).*couchstep + startcouch;
    case 'cradle'
        if Nshot>1
            error('Multi-shots scan is not supported in cradle mode.');
        end
        if ~isfield(SYS.protocol, 'CradleCurve')
            error('The simulation protocol of cradle scan shall include the field CradleCurve.');
        end
        [couch_z, turningpoints] = cradlecouch(SYS.protocol.CradleCurve, Nview, couchspeed, rotspeed, Nviewprot, startcouch);
    otherwise
        % static or topo
        couch_z = (0:Nview-1)'.*(couchspeed*inttime*1e-6) + ...
            (0:Nshot-1).*couchstep + startcouch;
end
couch = [zeros(Nview*Nshot, 1) -ones(Nview*Nshot, 1).*couchheight couch_z(:)];
% block
couch = couch(Startview : Startview + Nviewofcurrblk-1, :);

end


function [couch_z, turningpoints] = cradlecouch(CradleCurve, Nview, couchspeed, rotspeed, Nviewprot, startcouch)

if isfield(CradleCurve, 'curvedata') && isfile(CradleCurve.curvedata)
    % load reset curvedata
    curvedata = loaddata(CradleCurve.curvedata);
    couch_z = curvedata.couch_z(:) + startcouch;
    if length(couch_z) > Nview
        couch_z = couch_z(1:Nview);
    end
    turningpoints = curvedata.turningpoints;
    return
end

% reset cradle mode
switch CradleCurve.presetmode
    case 'braking'
        % braking mode
        t = ((0 : Nview-1).*(rotspeed/Nviewprot) - CradleCurve.brakingpoint) ./ CradleCurve.brakingtime;
        z = (-exp(2.*t)./2 + exp(2).*t + 1/2)./(exp(2)-1);
        z(t<=0) = t(t<=0);
        z(t>=1) = (exp(2)+1)/(exp(2)-1)/2;
        z = (z.*CradleCurve.brakingtime + CradleCurve.brakingpoint).*couchspeed;
        turningpoints = [floor(CradleCurve.brakingpoint/rotspeed*Nviewprot) + 1, ...
                         ceil((CradleCurve.brakingtime + CradleCurve.brakingpoint)/rotspeed*Nviewprot) + 1];
    case 'starting'
        % starting mode

        % TBC
    case 'twin'
        % starting-uniform-braking or starting-braking mode
        % TBC
    otherwise
        % error
end

couch_z = z(:) + startcouch;

end