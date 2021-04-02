% image sparse-net test

% head image1
imgpath = 'F:\data-Dier.Z\97_1\1.2.840.1.99.1.47.2.1609918020.1005';

[img0, dcminfo] = loaddicomimg(imgpath);
Nimg = size(img0, 3);

cut1 = 800;
img1 = img0.*(img0>cut1);

se = offsetstrel('ball', 5, 2);
img2 = zeros(size(img1));
for ii = 1:Nimg
    img2(:,:,ii) = imdilate(imerode(img1(:,:,ii), se), se);
end
% img2 = imdilate(

% % imLabel = bwlabel(imBw(:,:,ii));                %对各连通域进行标记
% stats = regionprops(imBw(:,:,ii),'Area');    %求各连通域的大小
% area = cat(1,stats.Area);   
