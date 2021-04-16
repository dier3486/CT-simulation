% image sparse-net test

% head image1
imgpath = 'F:\data-Dier.Z\TM\body\2.1608198353666.pd.2.1608199788691';

% load image
[img0, dcminfo] = loaddicomimg(imgpath);
Nimg = size(img0, 3);
hz = 2;

% imerode & imdilate
cut1 = 800;
img1 = img0.*(img0>cut1);
% se1 = offsetstrel('ball', 4, 2, 6);
% % se2 = offsetstrel('ball', 4, 2, 4);
% % se1 = strel('sphere',10);
% % se2 = strel('sphere',20);
% img2 = zeros(size(img1));
% for ii = 1:Nimg
%     img2(:,:,ii) = imdilate(imerode(img1(:,:,ii), se1), se1);
% %     img3(:,:,ii) = imerode(imdilate(imerode(img1(:,:,ii), se1), se2), se1);
% %     img3(:,:,ii) = imerode(imdilate(img2(:,:,ii), se1), se1);
% end
img2 = img1;

% region
cut2 = 800;
img3 = zeros(size(img1));
for ii = 1:Nimg
    imBw2 = img2(:,:,ii)>cut2;
    imLabel = bwlabel(imBw2);   
    stats = regionprops(imLabel,'Area');    %求各连通域的大小
    [~, imax] = max([stats.Area]);
    img3(:,:,ii) = imLabel==imax;
end

% blur
img4 = gaussblur(img3, 4).*img0;

Ybond = nan(2, 512, Nimg);
Ybond2 = Ybond;
cut4 = 800;
for ii = 1:Nimg
    for jj = 1:512
        Vjj = img4(:,jj,ii);
        index1 = find(Vjj>cut4, 1, 'first');
        if ~isempty(index1)
            Ybond(1, jj, ii) = index1 + (cut4-Vjj(index1))/(Vjj(index1)-Vjj(index1-1));
            index2 = find(Vjj>cut4, 1, 'last');
            Ybond(2, jj, ii) = index2 + (cut4-Vjj(index2))/(Vjj(index2+1)-Vjj(index2));
        end
    end
    Ybond2(1, :, ii) = smooth(Ybond(1, :, ii), 5);
    Ybond2(2, :, ii) = smooth(Ybond(2, :, ii), 5);
end
Xbond = zeros(2, Nimg);
Xbond(1,:) = smooth(sum(squeeze(isnan(Ybond(1,1:256,:))),1)+1);
Xbond(2,:) = smooth(512-sum(squeeze(isnan(Ybond(1,257:end,:))),1));
Xcut = 0.1;

Xbond2 = [1-Xcut Xcut; Xcut 1-Xcut]*Xbond;
% Nt1 = 256;
Nt1 = 64;
ht = 1/Nt1;
t1 = linspace(ht/2, 1-ht/2, Nt1)';
tgrid = [1-t1 t1]*Xbond2;
zgrid = repmat(1:Nimg, Nt1, 1);
Ysmp = zeros(Nt1, Nimg, 2);
Ysmp(:,:,1) = interp2(squeeze(Ybond2(1,:,:)), zgrid, tgrid);
Ysmp(:,:,2) = interp2(squeeze(Ybond2(2,:,:)), zgrid, tgrid);
Ysmp2 = Ysmp;
for ii = 1:Nt1
   Ysmp2(ii,:,1) = smooth(Ysmp(ii,:,1), 10);
   Ysmp2(ii,:,2) = smooth(Ysmp(ii,:,2), 10);
end

Ycent = mean(Ysmp2, 3);
for ii = 1:Nimg
    Ycent(:, ii) = smooth(Ycent(:, ii), 20);
end
for ii = 1:Nt1
    Ycent(ii, :) = smooth(Ycent(ii, :), 10);
end

nth = 16;
Ntheta = (2*nth+1)*2;

