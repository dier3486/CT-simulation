% nonlinear calibration test

% after cross talk correction
% fix center water again

% % load data
% datapath = 'E:\data\rawdata\bhtest\';
% datafile = {'df_w200c.mat', 'df_w300c.mat'};
% Nf = length(datafile);
% df = cell(1, Nf);
% for ifile = 1:Nf
%     df{ifile} = load(fullfile(datapath, datafile{ifile}));
% end
% 
% load non-linear1
corr_nl = load('E:\data\rawdata\bhtest\pnl_2.mat');
p0_nl = corr_nl.p_nl;
% 
% 
% Npixel = 864;
% Nslice = 64;
% Nview = 1152;

% % fix bh again
% m1 = 16;
% m2 = 48;
% cut = 1e4;
% 
% A1 = reshape(df{2}.rawdata, Npixel, []);
% B1 = reshape(df{2}.rawdata_bh, Npixel, []);
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
% 
% 
% 
% % non-linear2
% n = 2;
% p2_nl = zeros(Nslice, n);
% % islice = 5;
% t0 = zeros(1, n);
% t0(end) = 1;
% for islice = 1:Nslice
%     pindex = (1:Npixel) + (islice-1).*Npixel;
%     
%     x2 = iterpolyval(p0_nl(pindex, :), B2(pindex, :)./1000);
%     y2 = A2(pindex, :)./1000;
%     s2 = y2>0;
%     x2_s = double(x2(s2));
%     y2_s = double(y2(s2));
% 
% %     p1 = lsqnonlin(@(t) (iterpolyval(t, x2_s)-y2_s), [0 1]);
%     p2_nl(islice, :) = lsqnonlin(@(t) (iterinvpolyval(t, y2_s)-x2_s), t0);
%     t0 = p2_nl(islice, :);
% end
% p2_nl = repelem(p2_nl, Npixel, 1);
% save('E:\data\rawdata\bhtest\pnl2_1.mat', 'p2_nl');

% load non-linear2
corr_nl2 = load('E:\data\rawdata\bhtest\pnl2_1a.mat');
p2_nl = corr_nl2.p2_nl;

% fix nl again
m1 = 48;
m2 = 24;
cut = 1e4;

A1 = reshape(df{1}.rawdata, Npixel, []);
B1 = reshape(df{1}.rawdata_bh, Npixel, []);
A2 = reshape(df{2}.rawdata, Npixel, []);
B2 = reshape(df{2}.rawdata_bh, Npixel, []);
for iv = 1:Nview*Nslice
    ix_l = find(A1(:, iv)>cut/2, 1, 'first');
    ix_r = find(A1(:, iv)>cut/2, 1, 'last');
    A1([1:ix_l+m1-1 ix_r-m1+1:Npixel], iv) = 0;
    B1([1:ix_l+m1-1 ix_r-m1+1:Npixel], iv) = 0;
    ix_l = find(A2(:, iv)>cut, 1, 'first');
    ix_r = find(A2(:, iv)>cut, 1, 'last');
    A2([1:ix_l+m2-1 ix_r-m2+1:Npixel], iv) = 0;
    B2([1:ix_l+m2-1 ix_r-m2+1:Npixel], iv) = 0;
end
A1 = reshape(A1, [], Nview);
B1 = reshape(B1, [], Nview);
A2 = reshape(A2, [], Nview);
B2 = reshape(B2, [], Nview);

B1_fix = iterpolyval(p0_nl, B1./1000);
B1_fix = iterpolyval(p2_nl, B1_fix);

B2_fix = iterpolyval(p0_nl, B2./1000);
B2_fix = iterpolyval(p2_nl, B2_fix);

scl1 = reshape(A1./1000./B1_fix, Npixel, Nslice, Nview);
scl2 = reshape(A2./1000./B2_fix, Npixel, Nslice, Nview);

m_scl = 32;
Nmerge = 8;
Nslice_mg = Nslice/Nmerge;
alpha1 = 0.5;

p3_scl = ones(Npixel, Nslice_mg);
for ii = 1:Nslice_mg
    index_sl = (1:Nmerge) + (ii-1).*Nmerge;
    scl1_ii = squeeze(mean(scl1(:, index_sl, :), 2));
    scl2_ii = squeeze(mean(scl2(:, index_sl, :), 2));
    snan1 = sum(isnan(scl1_ii), 2)<=Nview/2;
    snan2 = sum(isnan(scl2_ii), 2)<=Nview/2;
    scl1_ii = mean(scl1_ii, 2, 'omitnan');
    scl2_ii = mean(scl2_ii, 2, 'omitnan');
    
    c1_l = find(snan1, 1, 'first');
    scl1_ii(1:c1_l-1) = mean(scl1_ii(c1_l:c1_l+m_scl-1));
    c1_r = find(snan1, 1, 'last');
    scl1_ii(c1_r+1:end) = mean(scl1_ii(c1_r-m_scl+1:c1_r));
    
    
    c2_l = find(snan2, 1, 'first');
    scl2_ii(1:c2_l-1) = mean(scl2_ii(c2_l:c2_l+m_scl-1));
    c2_r = find(snan2, 1, 'last');
    scl2_ii(c2_r+1:end) = mean(scl2_ii(c2_r-m_scl+1:c2_r));
    
    scl1_ii = scl1_ii-mean(scl1_ii)+mean(scl2_ii);
    p3_scl(:, ii) = scl1_ii.*alpha1+scl2_ii.*(1-alpha1);
end
p3_scl = repelem(p3_scl, 1, Nmerge);
p3_nl = p2_nl;
p3_nl(:,end) = p3_nl(:,end).*p3_scl(:);
save('E:\data\rawdata\bhtest\pnl3_1.mat', 'p3_nl');


