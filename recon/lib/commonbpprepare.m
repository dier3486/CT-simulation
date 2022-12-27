function recon = commonbpprepare(recon, protocol, system, BPprm)
% a subfunction of back projection prepare which will be used in
% backpojection or iteration recon prepare including the common prepare
% works of axial 2D/3D, helical and their interations
% prfmlow.recon = commonbpprepare(prfmlow.recon, prfmlow.protocol, prfmlow.system, BPprm);

if ~exist('BPprm', 'var')
    % while nargin<4
    BPprm = struct();
end

% FOV
if isfield(BPprm, 'FOV')
    recon.FOV = BPprm.FOV;
elseif isfield(protocol, 'reconFOV')
    recon.FOV = protocol.reconFOV;
else
    % default FOV = 500
    recon.FOV = 500;
end
% maxFOV
if isfield(BPprm, 'maxFOV')
    recon.maxFOV = BPprm.maxFOV;
elseif isfield(system, 'maxFOV')
    recon.maxFOV = system.maxFOV;
else
    % default maxFOV = 500 (or 650?)
    recon.maxFOV = 500.01;
end
% imagesize
if isfield(BPprm, 'imagesize')
    recon.imagesize = BPprm.imagesize;
elseif isfield(protocol, 'imagesize')
    recon.imagesize = protocol.imagesize;
else
    % default imagesize is 512x512
    recon.imagesize = [512, 512];
end
if length(recon.imagesize) == 1
    recon.imagesize = [recon.imagesize, recon.imagesize];
end
% image XY
if isfield(BPprm, 'center')
    recon.center = BPprm.center;
else
    recon.center = protocol.reconcenter;
end
% window
if isfield(BPprm, 'windowcenter')
    recon.windowcenter = BPprm.windowcenter;
else
    recon.windowcenter = protocol.windowcenter;
end
if isfield(BPprm, 'windowwidth')
    recon.windowwidth = BPprm.windowwidth;
else
    recon.windowwidth = protocol.windowwidth;
end
% kernel
if isfield(BPprm, 'kernel')
    recon.kernel = BPprm.kernel;
elseif isfield(protocol, 'reconkernel')
    recon.kernel = protocol.reconkernel;
else
    % default kernel?
    recon.kernel = '';
end
% imagethickness & imageincrement
if isfield(system, 'nominalslicethickness') && ~isempty(system.nominalslicethickness)
    % The real image/slice thickness could not exactly equal to the nominal thickness, 
    % e.g. when the nominal thickness could be 0.625 but the real thickness is 0.6222
    % and we truely move the couch by 0.6222*Nslice between each shot.
    zscale = system.detector.hz_ISO/system.nominalslicethickness;
    recon.imagethickness = protocol.imagethickness * zscale;
    recon.imageincrement = protocol.imageincrement * zscale;
    % mostly, we don't do that. If we do this that means we do the 0.6222 images undertable but show up 0.625 to customers.
    % Warning: like that idea?
else
    recon.imagethickness = protocol.imagethickness;
    recon.imageincrement = protocol.imageincrement;
end
% % tilt (moved to readrawdata.m)
% recon.gantrytilt = protocol.gantrytilt*(pi/180);
% couch (table) step & couch direction
recon.shotcouchstep = protocol.shotcouchstep;
recon.couchdirection = protocol.couchdirection;
% startcouch
recon.startcouch = protocol.startcouch;

% X-upsampling
if isfield(BPprm, 'upsampling') && ~isempty(BPprm.upsampling)
    recon.upsampling = BPprm.upsampling;
    % upsampling
    recon.delta_d_up = recon.delta_d/2;
    recon.midchannel_up = round(recon.midchannel*4-2)/2;
else
    recon.upsampling = false;
end
if isfield(BPprm, 'upsampgamma') && ~isempty(BPprm.upsampgamma)
    recon.upsampgamma = BPprm.upsampgamma;
else
    recon.upsampgamma = [0.7 0.8854];
end

% voxelsize
recon.voxelsize = single(recon.FOV/min(recon.imagesize));
% effective FOV
imageFOV = recon.FOV*sqrt(sum(recon.imagesize.^2))/min(recon.imagesize);
minFOV = 220;
recon.effFOV = max(min(imageFOV*0.85, recon.maxFOV), minFOV);

% XY grid of the image pixels
% xygrid = single((-(recon.imagesize-1)/2 : (recon.imagesize-1)/2).*recon.voxelsize);
% [X, Y] = ndgrid(xygrid);
xgrid = single((-(recon.imagesize(1)-1)/2 : (recon.imagesize(1)-1)/2).*recon.voxelsize);
ygrid = single((-(recon.imagesize(2)-1)/2 : (recon.imagesize(2)-1)/2).*recon.voxelsize);
[X, Y] = meshgrid(xgrid, ygrid);
Sxy = sqrt(X.^2+Y.^2) <= recon.effFOV/2;
% Sxy = sqrt(X.^2+Y.^2) <= inf;   % debug
recon.XY = ([X(Sxy) Y(Sxy)] - recon.center)./recon.SID;
recon.activeXY = Sxy;
recon.NactiveXY = size(recon.XY, 1);

% view angle and pitch (for helical)
switch lower(recon.scan)
    case 'axial'
        % It was
        %   1. in Axial 2D/3D recon viewangle = viewangle + startviewangle(ishot);
        %   2. in Axial interation viewangle = mod(viewangle(1:Nviewprot/2) + startviewangle(1), pi*2);
        recon.viewangle = (0:recon.Nviewprot-1)'.*recon.delta_view + recon.startviewangle(:)' + pi/2;
        recon.viewangle = single(recon.viewangle(:)');
    case 'helical'
        recon.viewangle = single((0:recon.Nview).*recon.delta_view + recon.startviewangle + pi/2);
        % pitch
        recon.pitchlength = abs(protocol.couchspeed*protocol.rotationspeed);
        recon.pitch = recon.pitchlength/(recon.Nslice*recon.delta_z);
        if isfield(protocol, 'pitchhelical')
            recon.nominalpitch = protocol.pitchhelical;
        else
            recon.nominalpitch = recon.pitchlength/(recon.Nslice*recon.imageincrement);
        end
    otherwise
        % static?
        recon.viewangle = recon.startviewangle + pi/2;
end
% rot 45
if isfield(BPprm, 'R45') && BPprm.R45
    % rot 45
    recon.viewangle = viewangle - pi/4;
end

if isfield(BPprm, 'viewblock') && ~isempty(BPprm.viewblock)
    recon.viewblock = BPprm.viewblock;
else
    recon.viewblock = recon.Nviewprot;
end

% Filter
if isfield(BPprm, 'Filter')
    % to put filter in recon
    if recon.upsampling
        recon.filter = loadfilter(BPprm.Filter, recon.Npixel*2, recon.delta_d/2);
    else
        recon.filter = loadfilter(BPprm.Filter, recon.Npixel, recon.delta_d);
    end
end

% Forward projection
recon.FPchannelpos = ((1:recon.Npixel)'-recon.midchannel).*recon.delta_d;
eta_C = recon.center(1).*sin(recon.viewangle) - recon.center(2).*cos(recon.viewangle);
indexstart_p = floor(recon.midchannel + (-recon.effFOV/2 + eta_C)./recon.delta_d);
% indexstart_p(indexstart_p<1) = 1;
% indexstart_n = ceil(midchannel*2 - indexstart_p);     % no
indexstart_n = floor(recon.midchannel*2 - indexstart_p);
% indexstart_n(indexstart_n>Npixel) = Npixel;
recon.indexstart = [indexstart_p(:)  indexstart_n(:)];
recon.effNp = ceil(recon.effFOV/recon.delta_d) + 2;

end