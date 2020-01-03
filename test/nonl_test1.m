% 

df1 = load('E:\data\rawdata\bhtest\df1.mat');   % small water
df2 = load('E:\data\rawdata\bhtest\df2.mat');   % big water
df3 = load('E:\data\rawdata\bhtest\df3.mat');   % no use

Npixel = 864;
Nslice = 64;
islice = 30;

A1 = squeeze(df1.rawdata(:,islice,:));
B1 = squeeze(df1.rawdata_bh(:,islice,:));    % small
A2 = squeeze(df2.rawdata(:,islice,:));
B2 = squeeze(df2.rawdata_bh(:,islice,:));    % big
A3 = squeeze(df3.rawdata(:,islice,:));
B3 = squeeze(df3.rawdata_bh(:,islice,:));

m = 16;
Nv = size(A1, 2);
cut = 1e4;


for ii = 1:Nv
    ix_l = find(B1(:,ii)>cut, 1, 'first');
    ix_r = find(B1(:,ii)>cut, 1, 'last');
    A1([1:ix_l+m-1 ix_r-m+1:Npixel], ii) = 0;
    B1([1:ix_l+m-1 ix_r-m+1:Npixel], ii) = 0;
    
    ix_l = find(B2(:,ii)>cut, 1, 'first');
    ix_r = find(B2(:,ii)>cut, 1, 'last');
    A2([1:ix_l+m-1 ix_r-m+1:Npixel], ii) = 0;
    B2([1:ix_l+m-1 ix_r-m+1:Npixel], ii) = 0;
end

n = 2;
p = zeros(Npixel, n);
p1 = zeros(N2-N1, n);
Nt = 100;
d = 2;
w_cut = 5;
w2 = 0.3;   % weight for big water
options = optimoptions('lsqnonlin','Display','off');
% for ii = N1:N2
for ii = 1:Npixel
    x = double([B1(ii, :) B2(ii, :)]./1000);
    y = double([A1(ii, :) A2(ii, :)]./1000);
    if sum(y>0)<100
        continue;
    end
    miny = min(y(y>0));
    miny = max(5, miny-10);
    maxy = max(y) + 5;
    ty = linspace(miny, maxy, Nt);
    ty = ty(:);
    winf = exp(-(ty-y).^2./d^2);
    w = sum(winf,2);
    tx = (winf*x(:))./w;
    swa = ty>max(A1(ii, :))/1000;
    w(swa) = w(swa).*w2;
    sw = w>w_cut;
    tx = tx(sw);
    ty = ty(sw);
    w = w(sw);
    pii = lsqnonlin(@(t) (polyval(t, tx).*tx-ty).*sqrt(w), [1 1], [], [], options);
    p(ii, :) = pii;
end




% for ii = 400:400
%     x = B1(ii, :);
%     y = A1(ii, :);
%     sy = y>0;
%     if sum(sy)<40
%         continue;
%     end
%     x1 = x(sy);
%     y1 = y(sy);
%     [x1, sx] = sort(x1);
%     y1 = y1(sx);
%     
%     p(ii, :) = polyfit(x1, y1./x1, n-1);
% end