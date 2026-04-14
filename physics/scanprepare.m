function [samplekeV, viewangle, couch, shotindex, gantrytilt, couchgear] = scanprepare(SYS, iblock)
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
    case {'axial', 'helical', 'conveyor'}
        % rotation
        viewangle = (0 : round(Nrot*Nviewprot-1)).*(pi*2 / Nviewprot).*rotdirect;
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
        couchgear = zeros(1, Nview*Nshot);
    case 'helical'
        couch_z = (0:Nview-1)'.*(couchspeed*rotspeed/Nviewprot) + ...
            (0:Nshot-1).*couchstep + startcouch;
        couchgear = ones(1, Nview*Nshot).*sign(couchspeed);
    case 'conveyor'
        if Nshot>1
            error('Multi-shots scan is not supported in conveyor mode.');
        end
        if isfield(SYS.protocol, 'ConveyorCurve')
            ConveyorCurve = SYS.protocol.ConveyorCurve;
        elseif isfield(SYS.protocol, 'CradleCurve')
            ConveyorCurve = SYS.protocol.CradleCurve;
        else
            error('The simulation protocol of conveyor scan shall include the field ConveyorCurve.');
        end
        [couch_z, couchgear] = conveyorcouch(ConveyorCurve, Nview, couchspeed, rotspeed, Nviewprot, startcouch);
    otherwise
        % static or topo
        couch_z = (0:Nview-1)'.*(couchspeed*inttime*1e-6) + ...
            (0:Nshot-1).*couchstep + startcouch;
end
couch = [zeros(Nview*Nshot, 1) -ones(Nview*Nshot, 1).*couchheight couch_z(:)];
% block
couch = couch(Startview : Startview + Nviewofcurrblk-1, :);
couchgear = couchgear(Startview : Startview + Nviewofcurrblk-1);

end

function [couch_z, couchgear] = conveyorcouch(ConveyorCurve, Nview, couchspeed, rotspeed, Nviewprot, startcouch)

if isfield(ConveyorCurve, 'curvedata') && isfile(ConveyorCurve.curvedata)
    % load reset curvedata
    curvedata = loaddata(ConveyorCurve.curvedata);
    couch_z = curvedata.couch_z(:) + startcouch;
    if length(couch_z) > Nview
        couch_z = couch_z(1:Nview);
    end
    couchgear = turning2gear(Nview, couchspeed, ConveyorCurve.turningpoints, ConveyorCurve.presetmode);
    return
end

% reset conveyor mode
switch ConveyorCurve.presetmode
    case 'braking'
        % braking mode -- starting form uniform speed
        couchgear = ones(1, Nview);
        t = (0 : Nview-1).*(rotspeed/Nviewprot);
        z = t.*couchspeed;
        brakingflag = 1;
    case 'starting'
        % starting mode -- starting form static
        couchgear = zeros(1, Nview);
        t = (0 : Nview-1).*(rotspeed/Nviewprot);
        z = t.*0;
        brakingflag = 0;
    otherwise
        error('Unkown presetmode %s.', ConveyorCurve.presetmode);
end

brakingpoint = ConveyorCurve.brakingpoint(:);
brakingtime = ConveyorCurve.brakingtime(:);
Nb = length(brakingpoint);
if length(brakingtime) < Nb
    brakingtime = repmat(brakingtime, Nb/length(brakingtime), 1);
end

z0 = 0; t0 = 0;
for ii = 1:Nb
    if ii > 1
        t_dag = (t0 - brakingpoint(ii))/brakingtime(ii-1);
        t_dag = max(t_dag, 0);
    else
        t_dag = 0;
    end
    t_ii = (t - brakingpoint(ii))./brakingtime(ii) + t_dag;
    s2 = t_ii>=1;
    s1 = t_ii > t_dag & ~s2;

    if mod(ii, 2) == brakingflag
        % braking
        z_dag = (-exp(2*t_dag)./2 + exp(2)*t_dag + 1/2)./(exp(2)-1);
        z1 = (-exp(2.*t_ii(s1))./2 + exp(2).*t_ii(s1) + 1/2)./(exp(2)-1) - z_dag*2;
        z2 = (exp(2)+1)/(exp(2)-1)/2 - z_dag*2;

        z_eps = eps*2/brakingtime(ii);
        if ~isempty(z1) && (z2 - z1(end) < z_eps)
            z1(end) = z1(end) - z_eps;
        end

        z(s1) = (z1.*brakingtime(ii) + z0 + max(brakingpoint(ii)-t0, 0)).*couchspeed;

        z0 = z0 + max(brakingpoint(ii)-t0, 0) + z2*brakingtime(ii);
        t0 = brakingpoint(ii) + brakingtime(ii) - t_dag*brakingtime(ii);

        z(s2) = z0.*couchspeed;

        couchgear(s1) = 3;      % slowing down
        couchgear(s2) = 0;      % stopped

    else
        % starting
        z_dag = (exp(-2*t_dag)/2 + t_dag - 1/2)/(1-exp(-2));
        z1 = (exp(-2.*t_ii(s1))./2 + t_ii(s1) - 1/2)./(1-exp(-2)) - z_dag*2;
        z2 = (exp(2)+1)/(exp(2)-1)/2 - z_dag*2;

        z_eps = eps*2/brakingtime(ii);
        if ~isempty(z1) && (z1(1) < z_eps)
            z1(1) = z1(1) + z_eps;
        end

        z(s1) = (z1.*brakingtime(ii) + z0).*couchspeed;

        z0 = z0 + z2*brakingtime(ii);
        t0 = brakingpoint(ii) + brakingtime(ii) - t_dag*brakingtime(ii);

        z(s2) = (z0 + t(s2) - t0).*couchspeed;

        couchgear(s1) = 2;      % speeding up
        couchgear(s2) = 1;      % speeded
    end
        
end

couchgear = couchgear.*sign(couchspeed);
couch_z = z(:) + startcouch;

end


function couchgear = turning2gear(Nview, couchspeed, turningpoints, presetmode)

couchgear = ones(1, Nview).*sign(couchspeed);

switch lower(presetmode)
    case 'braking'
        % braking mode
        for ii = 1:size(turningpoints, 2)
            if mod(ii, 2)
                couchgear(turningpoints(1, ii)+1:turningpoints(2, ii)-1) = sign(couchspeed)*3;
                couchgear(turningpoints(2, ii):end) = 0;
            else
                couchgear(turningpoints(1, ii)+1:turningpoints(2, ii)-1) = sign(couchspeed)*2;
                couchgear(turningpoints(2, ii):end) = sign(couchspeed);
            end
        end
    case 'starting'
        % starting mode
        % TBC
    otherwise
        % error
end
end
