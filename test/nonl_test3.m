% nonlinear calibration test

% after cross talk correction
% fix big water again

% % load data
% datapath = 'E:\data\rawdata\bhtest\';
% datafile = 'df2_w300c.mat';
% df = load(fullfile(datapath, datafile));

% Npixel = 864;
% Nslice = 64;
% Nview = 1152;
% 
% m1 = 16;
% m2 = 48;
% cut = 1e4;
% 
% A1 = reshape(df.rawdata, Npixel, []);
% B1 = reshape(df.rawdata_bh, Npixel, []);
% A2 = A1;  B2 = B1;
% for iv = 1:Nview*Nslice
%     ix_l = find(A1(:, iv)>cut, 1, 'first');
%     ix_r = find(A1(:, iv)>cut, 1, 'last');
%     A1([1:ix_l+m1-1 ix_r-m1+1:Npixel], iv) = 0;
%     B1([1:ix_l+m1-1 ix_r-m1+1:Npixel], iv) = 0;
%     A2([1:ix_l+m2-1 ix_r-m2+1:Npixel], iv) = 0;
%     B2([1:ix_l+m2-1 ix_r-m2+1:Npixel], iv) = 0;
% end
% A1 = reshape(A1, [], Nview);
% B1 = reshape(B1, [], Nview);
% A2 = reshape(A2, [], Nview);
% B2 = reshape(B2, [], Nview);

% % load non-linear1
% corr_nl = load('E:\data\rawdata\bhtest\pnl_2.mat');
% p0_nl = corr_nl.p_nl;

% non-linear2
n = 2;
p2_nl = zeros(Nslice, n);
% islice = 5;
t0 = zeros(1, n);
t0(end) = 1;
for islice = 1:Nslice
    pindex = (1:Npixel) + (islice-1).*Npixel;
    
    x2 = iterpolyval(p0_nl(pindex, :), B2(pindex, :)./1000);
    y2 = A2(pindex, :)./1000;
    s2 = y2>0;
    x2_s = double(x2(s2));
    y2_s = double(y2(s2));

%     p1 = lsqnonlin(@(t) (iterpolyval(t, x2_s)-y2_s), [0 1]);
    p2_nl(islice, :) = lsqnonlin(@(t) (iterinvpolyval(t, y2_s)-x2_s), t0);
    t0 = p2_nl(islice, :);
end
p2_nl = repelem(p2_nl, Npixel, 1);


% x2_fix1 = iterpolyval(p1, x2);
% x2_fix2 = iterpolyval(p2, x2);
% c2 = mean(y2./x2, 2, 'omitnan');
% c2_fix1 = mean(y2./x2_fix1, 2, 'omitnan');
% c2_fix2 = mean(y2./x2_fix2, 2, 'omitnan');
% c2(sum(s2,2)<Nview/2) = nan;
% c2_fix1(sum(s2,2)<Nview/2) = nan;
% c2_fix2(sum(s2,2)<Nview/2) = nan;
% 
% figure;
% plot([c2_fix1 c2_fix2]);
% 
% yy = linspace(0, max(y2(:)), 500);
% xx_p1 = iterinvpolyval(p1, yy);
% xx_p2 = iterinvpolyval(p2, yy);
% figure; hold on
% plot(x2(:).*K, x2(:)./y2(:), '.');
% plot(xx_p1.*K, xx_p1./yy);
% plot(xx_p2.*K, xx_p2./yy);


