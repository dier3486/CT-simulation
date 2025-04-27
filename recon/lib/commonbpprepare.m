function recon = commonbpprepare(recon, protocol, system, BPprm)
% a subfunction of back projection prepare which will be used in
% backpojection or iteration recon prepare including the common prepare
% works of axial 2D/3D, helical and their interations
% prfmlow.recon = commonbpprepare(prfmlow.recon, prfmlow.protocol, prfmlow.system, BPprm);

if ~exist('BPprm', 'var')
    % while nargin<4
    BPprm = struct();
end

% scan
recon.scan = protocol.scan;

% Nview
% I know the Nview was prepared by rebin
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
% effFOV
if isfield(BPprm, 'effectiveFOV')
    recon.effFOV = min(BPprm.effectiveFOV, recon.maxFOV);
elseif isfield(protocol, 'effectiveFOV')
    recon.effFOV = min(protocol.effectiveFOV, recon.maxFOV);
end
% imagesize
if isfield(BPprm, 'imagesize')
    recon.imagesize = BPprm.imagesize(:)';
elseif isfield(protocol, 'imagesize')
    recon.imagesize = protocol.imagesize(:)';
else
    % default imagesize is 512x512
    recon.imagesize = [512, 512];
end
if length(recon.imagesize) == 1
    recon.imagesize = [recon.imagesize, recon.imagesize];
end
% recon center
if isfield(BPprm, 'center')
    recon.center = BPprm.center(:)';
else
    recon.center = protocol.reconcenter(:)';
end
% the recon center is the ISO center postion on image, therefore the position of image center to th ISO is -recon.center.
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
    % re-scale the imagethickness and imageincrement in protocol
    zscale = system.detector.hz_ISO/system.nominalslicethickness;
    recon.imagethickness = protocol.imagethickness * zscale;
    recon.imageincrement = protocol.imageincrement * zscale;
    % The real image/slice thickness could not exactly equal to the nominal thickness, 
    % e.g. the nominal thickness could be 0.625 but the real thickness is 0.6222
    % and we truely move the couch by 0.6222*Nslice between each shot. 
    % The '0.625' is used to fool the customers and the 0.6222 is what really happened.
    % mostly, we don't do that. If we do that it means we do the 0.6222 images undertable but show up 0.625 to customers.
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
% image instance number start
recon.InstanceStart = 0;

% upsampled?
if ~isfield(recon, 'upsampled')
    recon.upsampled = false;
end
% X-upsampling
if ~recon.upsampled
    if isfield(BPprm, 'upsampling') && ~isempty(BPprm.upsampling)
        recon.upsampling = BPprm.upsampling;
    else
        % default upsampling flag is true while the BP data has not been filterd, or it is false.
        recon.upsampling = ~recon.filtered;
    end
    if recon.upsampling
        % upsampling
        recon.Npixel_up = recon.Npixel*2;
        recon.delta_d_up = recon.delta_d/2;
        recon.midchannel_up = round(recon.midchannel*4-2)/2;
    else
        recon.Npixel_up = recon.Npixel;
        recon.delta_d_up = recon.delta_d;
        recon.midchannel_up = recon.midchannel;
    end
    if isfield(BPprm, 'upsampgamma') && ~isempty(BPprm.upsampgamma)
        recon.upsampgamma = BPprm.upsampgamma;
    else
        recon.upsampgamma = [0.7 0.8854];
    end
else
    % has been upsampled
    recon.upsampling = false;
end


% voxelsize
if isfield(BPprm, 'voxelsize') && ~isempty(BPprm.voxelsize)
    recon.voxelsize = BPprm.voxelsize;
else
    recon.voxelsize = single(recon.FOV/min(recon.imagesize));
end
% effective FOV
if ~isfield(recon, 'effFOV')
    imageFOV = recon.FOV*sqrt(sum(recon.imagesize.^2))/min(recon.imagesize);
    minFOV = 220;
    recon.effFOV = max(min(imageFOV*0.85, recon.maxFOV), minFOV);
end

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

% pitch (for helical)
if any(strcmpi(recon.scan, {'helical', 'cradle'}))
    recon.pitchlength = abs(protocol.couchspeed*protocol.rotationspeed);
    recon.pitch = recon.pitchlength/(recon.Nslice*recon.delta_z);
    if isfield(protocol, 'pitchhelical')
        recon.nominalpitch = protocol.pitchhelical;
    else
        recon.nominalpitch = recon.pitchlength/(recon.Nslice*recon.imageincrement);
    end
end

% rot 45
if isfield(BPprm, 'R45') && BPprm.R45
    % rot 45
    recon.Rot45 = true;
else
    recon.Rot45 = false;
end

if isfield(BPprm, 'viewblock') && ~isempty(BPprm.viewblock)
    recon.viewblock = BPprm.viewblock;
else
    recon.viewblock = ceil(recon.Nviewprot / 2);
end

% Filter
if isfield(BPprm, 'Filter')
    % to put filter in recon
    recon.filter.basicfilter = loadfilter(BPprm.Filter, recon.Npixel_up, recon.delta_d_up);
    % 'alpha' mixing (mostly for iteratin)
    if isfield(BPprm.Filter, 'alpha')
        % baseline ram-lak
        filterRL = filterdesign('ram-lak', recon.Npixel_up, recon.delta_d_up, inf);
        % mix with ram-lak
        recon.filter.basicfilter = single(recon.filter.basicfilter.*(1-BPprm.Filter.alpha) + filterRL.*BPprm.Filter.alpha);
    end
elseif ~isfield(recon, 'filtered') || ~recon.filtered
    % load a default filter
    defaultFilter.name = 'none';
    recon.filter.basicfilter = loadfilter(defaultFilter, recon.Npixel_up, recon.delta_d_up);
end

% is concyclic (detector frame)
if isfield(BPprm, 'concyclic')
    recon.concyclic = BPprm.concyclic;
elseif isfield(system.detector, 'concyclic')
    recon.concyclic = system.detector.concyclic;
else
    recon.concyclic = false;
end

% Forward projection
recon.forward = struct();
recon.forward.FPchannelpos = ((1:recon.Npixel)'-recon.midchannel).*recon.delta_d;
recon.forward.effNp = ceil(recon.effFOV/recon.delta_d) + 2;

% viewread/imagewritten
recon.viewread = 0;
recon.imagewritten = 0;
recon.availwritten = 0;

end