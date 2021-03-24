function [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status)
% recon node, Axial 2D BP (afrer reconnode_BPprepare)
% [dataflow, prmflow, status] = reconnode_Axial2DBackprojection(dataflow, prmflow, status);

% BP parameters
BPprm = prmflow.pipe.(status.nodename);
% interp method (linear for 4points)
if isfield(BPprm, 'interp')
    interpmethod = BPprm.interp;
else
    interpmethod = '4points';
end

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end

% recon parameters
SID = prmflow.recon.SID;
Nshot = prmflow.recon.Nshot;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
startviewangle = prmflow.recon.startviewangle + pi/2;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d/SID;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Neighb = prmflow.recon.Neighb;
Nextslice = prmflow.recon.Nextslice;
Nimage = prmflow.recon.Nimage;
reconcenter = prmflow.recon.center;
% imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
blocksize = 128^2;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);

% ini image
dataflow.image = zeros(imagesize, imagesize, Nimage, 'single');

% XY
h = FOV/imagesize/SID;
xygrid = single((-(imagesize-1)/2 : (imagesize-1)/2).*h);
[X, Y] = ndgrid(xygrid);
XY = [X(:) Y(:)] - reconcenter./SID;
% view angle
viewangle_prot = single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot));

% Chninterp  
Chninterp.delta_d = delta_d;
Chninterp.midchannel = midchannel;
Chninterp.channelindex = single(1:Npixel)';

% Zinterp (it was prepared in reconnode_BPprepare) 
Zinterp = prmflow.recon.Zinterp;

% to GPU
if GPUonoff
    XY = gpuArray(XY);
    viewangle_prot = gpuArray(viewangle_prot);
    Chninterp = everything2single(Chninterp, 'any', 'gpuArray');
    Zinterp = everything2single(Zinterp, 'any', 'gpuArray');
end

for ishot = 1:Nshot
    % image index and shotflag
    [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot, couchdirection);
    % view angles
    viewangle = viewangle_prot + startviewangle(ishot);  
    % 3D BP
    img_shot = backproj3Dslope(dataflow.rawdata(:, :, :, ishot), viewangle, XY, Neighb, Chninterp, Zinterp, shotflag, ...
        interpmethod, blocksize, GPUonoff);
    % get img
    dataflow.image(:,:,imageindex) = dataflow.image(:,:,imageindex) + reshape(gather(img_shot), imagesize, imagesize, []);
end
% normalize by Nviewprot, (pi/2 due to the filter)
dataflow.image = dataflow.image.*(pi/Nviewprot/2);

% reorderflag = couchdirection < 0;
% dataflow.image = imagereorder(dataflow.image, Nslice, reorderflag);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [imageindex, shotflag] = getimageindex(Nslice, Nextslice, Nshot, ishot, couchdirection)

if couchdirection<0
    % backward couch
    imageindex = (1:Nextslice) + ((ishot-1)*Nslice-(Nextslice-Nslice)/2);
    if Nshot==1
        shotflag = 3;
    elseif ishot==1
        imageindex(1:(Nextslice-Nslice)/2) = [];
        shotflag = 1;
    elseif ishot==Nshot
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        shotflag = 2;
    else
        shotflag = 0;
    end
else
    % forward couch
    imageindex = (1:Nextslice) + ((Nshot-ishot)*Nslice-(Nextslice-Nslice)/2);
    if Nshot==1
        shotflag = 3;
    elseif ishot==Nshot
        imageindex(1:(Nextslice-Nslice)/2) = [];
        shotflag = 1;
    elseif ishot==1
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        shotflag = 2;
    else
        shotflag = 0;
    end
end

end