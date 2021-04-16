% image sparse-net test

% head image1
imgpath = 'F:\data-Dier.Z\97_1\1.2.840.1.99.1.47.2.1609918020.1005';

% load image
[img0, dcminfo] = loaddicomimg(imgpath);
Nimg = size(img0, 3);
hz = 2;

% imerode & imdilate
cut1 = 800;
img1 = img0.*(img0>cut1);
se1 = offsetstrel('ball', 4, 2, 6);
% se2 = offsetstrel('ball', 4, 2, 4);
% se1 = strel('sphere',10);
% se2 = strel('sphere',20);
img2 = zeros(size(img1));
for ii = 1:Nimg
    img2(:,:,ii) = imdilate(imerode(img1(:,:,ii), se1), se1);
%     img3(:,:,ii) = imerode(imdilate(imerode(img1(:,:,ii), se1), se2), se1);
%     img3(:,:,ii) = imerode(imdilate(img2(:,:,ii), se1), se1);
end

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

% r-theta
Nth = 128;
rawrth = rthetatrans(img3, [0 0], Nth, 1.0);
rawrth = fillmissing(rawrth, 'constant', 0);

% bnd
Nraw = size(rawrth, 1);
crth = cumsum(rawrth);
Cbnd = [Nraw/2-squeeze(sum(crth==crth(end,:,:))); Nraw/2-squeeze(sum(crth==0))];

% smooth xy
span_th = 0.05;
edge = 32;
Cbnd2 = Cbnd;
for ii = 1:Nimg
    tmp = [Cbnd(end-edge+1:end, ii); Cbnd(:,ii); Cbnd(1:edge, ii)];
    tmp = smooth(tmp, span_th, 'loess');
    Cbnd2(:, ii) = tmp(edge+1:end-edge);
end
% smooth z
span_z = 0.03;
for ii = 1:Nth*2
    tmp = [repmat(Cbnd2(ii, 1), 1, edge) Cbnd2(ii,:) repmat(Cbnd2(ii, end), 1, edge)];
    tmp = smooth(tmp, span_z, 'loess');
    Cbnd2(ii, :) = tmp(edge+1:end-edge);
end

% sequence
seqset = sobolset(3, 'skip', 0);
Nseq = 10000; 
Pseq = net(seqset, Nseq);
% [u, v] = square2circle(Pseq(:,1).*2-1, Pseq(:,2).*2-1, true);
v = samplesetinobject(Pseq, 'spcylinder', true);

% z
Sc = sum(Cbnd2.^2).*(pi/Nth/2);
Scum = [0 cumsum(Sc)./sum(Sc)];
Zsample = (0:Nimg).*hz;
z = interp1(Scum, Zsample, (v(:,3)+1)./2);

% theta
Bndcum = cumsum(Cbnd2, 1);
Bndcum = [zeros(1, Nimg); Bndcum./Bndcum(end, :);];
thetasample = linspace(0, pi*2, Nth*2+1);
Bndsample = linspace(0, 1, Nth*2+1);
Bndcurve = zeros(Nth*2+1, Nimg+2);
for ii = 1:Nimg
    Bndcurve(:, ii+1) = interp1(Bndcum(:, ii), thetasample, Bndsample);
end
Bndcurve(:, 1) = Bndcurve(:, 2);
Bndcurve(:, end) = Bndcurve(:, end-1);
Zgrid = (-1/2:Nimg+1/2).*hz;
theta = interp2(Zgrid, Bndsample, Bndcurve, z, v(:,2)./(pi*2));

% r
Cbndsample = [Cbnd2; Cbnd2(1, :)];
Cbndsample = [Cbndsample(:, 1) Cbndsample Cbndsample(:, end)];
r = v(:,1).*interp2(Zgrid, thetasample, Cbndsample, z, theta);

% x y
x = r.*cos(theta);
y = r.*sin(theta);

% plot
thc = (0:Nth*2-1)'.*(pi/Nth);
Xc = Cbnd2.*cos(thc);    Xc = [Xc; Xc(1,:)];
Yc = Cbnd2.*sin(thc);    Yc = [Yc; Yc(1,:)];
Zc = repmat(((1:Nimg)-0.5).*hz, Nth*2+1, 1);
figure;
hold on;
mesh(Xc, Yc, Zc);
plot3(x, y, z, 'r.');
axis equal
hidden off
