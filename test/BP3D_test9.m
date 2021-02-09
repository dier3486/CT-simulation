% BP test code
% for 3D Axial slope rebin, shots fill up
% final test

% load('E:\data\simulation\TM\test\tilt_test3.mat');

% inputs are dataflow, prmflow
if exist('df0', 'var')
    dataflow = df0;
    clear df0;
end

if exist('pf0', 'var')
    prmflow = pf0;
    clear pf0;
end

% detector = prmflow.system.detector;
% SID = detector.SID;
% SDD = detector.SDD;
% delta_z = detector.hz_ISO;
SID = prmflow.recon.SID;
delta_z = prmflow.recon.delta_z;
Nshot = prmflow.recon.Nshot;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
startviewangle = prmflow.recon.startviewangle + pi/2;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d/SID;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nimage = prmflow.recon.Nimage;
% imageincrement = prmflow.recon.imageincrement;
% imageincrement = delta_z;
imageincrement = 0.6;   % imageincrement<delta_z
reconcenter = prmflow.recon.center;
% imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
gantrytilt = single(prmflow.recon.gantrytilt);

FOV = 300;

gpuDevice;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nviewprot, Nshot);

% ini image
img = zeros(imagesize, imagesize, Nimage, 'single');

xygrid = single((-(imagesize-1)/2 : (imagesize-1)/2).*(h/SID));
[X, Y] = ndgrid(xygrid);
XY = gpuArray([X(:) Y(:)] - reconcenter);

viewangle_prot = gpuArray(single(linspace(0, pi*2-pi*2/Nviewprot, Nviewprot)));

% Chninterp
Chninterp.delta_d = delta_d;
Chninterp.midchannel = midchannel;
Chninterp.channelindex = single(1:Npixel)';
% to GPU
Chninterp = everything2single(Chninterp, 'any', 'gpuArray');

% Zinterp
gamma_coeff = [0.6, 1.4];
Nzsample = [512 512];
maxFOV = 500;
Zinterp = omiga4table(gamma_coeff, Nzsample, maxFOV, FOV, SID, Nslice, gantrytilt);
% to GPU
Zinterp = everything2single(Zinterp, 'any', 'gpuArray');

% Nshot = 1;
for ishot = gpuArray(1:Nshot)
% for ishot = 1:1
tic;
    imageindex = (1:Nextslice) + ((ishot-1)*Nslice-(Nextslice-Nslice)/2);
    gatherindex = 1:Nextslice;
    shotflag = 0;
    if ishot==1
        imageindex(1:(Nextslice-Nslice)/2) = [];
        gatherindex(1:(Nextslice-Nslice)/2) = [];
        shotflag = 1;
    end
    if ishot==Nshot
        imageindex(end-(Nextslice-Nslice)/2+1:end) = [];
        gatherindex(end-(Nextslice-Nslice)/2+1:end) = [];
        shotflag = 2;
    end
    if Nshot==1
        shotflag = 3;
    end

    viewangle = viewangle_prot + startviewangle(ishot);  
    
    img_shot = backproj3Dslope(dataflow.rawdata(:, :, :, ishot), viewangle, XY, Nedge, Chninterp, Zinterp, shotflag, '4points', true);
    % get img
    img(:,:,imageindex) = img(:,:,imageindex) + reshape(gather(img_shot), imagesize, imagesize, []);
toc;
end

img = img.*(pi/Nviewprot/2);

