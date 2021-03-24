inputpath = 'F:\data-Dier.Z\TL\headBC\img\202';
infiles = dir(fullfile(inputpath, '*.dcm'));

Nimg = size(infiles, 1);
img = cell(Nimg, 1);

for ii = 1:Nimg
    dcm(ii) = dicominfo(fullfile(inputpath, infiles(ii).name));
    img{ii}= dicomread(fullfile(inputpath, infiles(ii).name));
end

h= 0.625;
gantrytilt = -6.5;
outputpath = 'F:\data-Dier.Z\TL\headBC\img\202c';
dcm1 = dcm;
img1 = img;
thetatilt = gantrytilt/180*pi;
% A1=[1 0 0 0 1 0]';
A1 = [1 0 0 0 cos(thetatilt) sin(thetatilt)]';
% A1 = [0 cos(thetatilt) sin(thetatilt) 1 0 0 ]';
% A1=[1 0 0 0 1 0]';
for ii = 1:Nimg
%     dcm1(ii).ImageOrientationPatient = A1;
%     dcm1(ii).GantryDetectorTilt = gantrytilt/180*pi;
%     dcm1(ii).ImagePositionPatient(3) = ii*h*cos(thetatilt);
%     dcm1(ii).SliceLocation = -(mod(ii,32)-16.5)*(h*cos(thetatilt));
%     dcm1(ii).SliceLocation = 0;
%     dcm1(ii).PixelSpacing = dcm1(ii).PixelSpacing.*[cos(thetatilt); 1];
%     dcm1(ii).PixelAspectRatio = [cos(thetatilt); 1];
%     dcm1(ii).SliceThickness = h*cos(thetatilt);
%     img1{ii}(:) = 1;
%     img1{ii}(250:280,150:400) = 4000;
    dcm1(ii).RescaleSlope = 0.5;
    dcm1(ii).RescaleIntercept = -512;
    dicomwrite(img1{ii}, fullfile(outputpath, infiles(ii).name), dcm1(ii), 'CreateMode','copy');
end