% BP test code
% for 3D Axial no tilt, shots fill up

% inputs are dataflow, prmflow
if exist('df0', 'var')
    dataflow = df0;
    clear df0;
end

if exist('pf0', 'var')
    prmflow = pf0;
    clear pf0;
end

detector = prmflow.system.detector;
SID = detector.SID;
SDD = detector.SDD;
delta_z = detector.hz_ISO;
Nshot = prmflow.recon.Nshot;
FOV = prmflow.recon.FOV;
Nviewprot = prmflow.recon.Nviewprot;
startviewangle = prmflow.recon.startviewangle;
imagesize = prmflow.recon.imagesize;
midchannel = prmflow.recon.midchannel;
delta_d = prmflow.recon.delta_d;
Nslice = prmflow.recon.Nslice;
Nimage = prmflow.recon.Nimage;
% imageincrement = prmflow.recon.imageincrement;
imageincrement = delta_z;
imagecenter = prmflow.recon.imagecenter;
couchdirection = prmflow.recon.couchdirection;
% set center to (0, 0)
imagecenter(:, 1:2) = 0;
% no tilt

h = FOV/imagesize;


% ini image
image = zeros(imagesize, imagesize, Nimage);

xygrid = (-(imagesize-1)/2 : (imagesize-1)/2).*h;
[X, Y] = ndgrid(xygrid);

zgrid = (-(Nslice-1)/2 : (Nslice-1)/2);
% slicegrid = (-(Nslice-1)/2 : (Nslice-1)/2).*delta_z;
midslice = (Nslice+1)/2;

% for ishot = 1:Nshot
for ishot = 1:1
    sliceindex = (1:Nslice) + (ishot-1)*Nslice;
    viewangle = startviewangle(ishot) + linspace(0, pi*2, Nviewprot);
%     viewangle = linspace(0, pi*2, Nviewprot);
    costheta = cos(viewangle);
    sintheta = sin(viewangle);
%     Xis = X(:) - imagecenter(sliceindex, 1)';
%     Yis = Y(:) - imagecenter(sliceindex, 2)';
    for iview = 1:Nviewprot
        Eta = -X.*sintheta(iview) + Y.*costheta(iview);
        Zeta = X.*costheta(iview) + Y.*sintheta(iview);
        t_chn = Eta./delta_d + midchannel;
        
        D = sqrt(SID^2 - Eta.^2);
        t1 = (D + Zeta).*(delta_z/imageincrement/SID);
        t2 = (D - Zeta).*(delta_z/imageincrement/SID);
        gap = Nslice - (t1+t2).*(Nslice-1)./2;
        
        kg1 = (reshape(1:Nslice/2,1,1,[])-1/2)./t1+1/2;
        kg2 = (reshape(1:Nslice/2,1,1,[])-1/2-Nslice)./t2+1/2 + Nslice;
        
        s1 = kg1<Nslice/2;
        s2 = (kg2>Nslice/2+1) & ~s1;
        s_gap = ~s1 & ~s2;
        
        h_z = zeros(size(kg1));
        t_z = zeros(size(kg1));
        t_z(s1) = kg1(s1);
        t_z(s2) = kg2(s2);
        
        t1 = repmat(t1, 1, 1, Nslice/2);
        t2 = repmat(t2, 1, 1, Nslice/2);
        gap = repmat(gap, 1, 1, Nslice/2);
        t_z(s_gap) = Nslice/2 + (kg1(s_gap) - Nslice/2).*t1(s_gap)./gap(s_gap);
        
        h_z(s1) = t1(s1);
        h_z(s2) = t2(s2);
        h_z(s_gap) = gap(s_gap);
        
        
    end
end

