% nonlinear calibration test

% after cross talk correction
% fix center water again

% % load data
% datapath = 'E:\data\rawdata\bhtest\';
% datafile = {'df2_w200c.mat', 'df2_w200off100.mat', 'df2_w300c.mat', 'df2_w300off100.mat'};
% Nf = length(datafile);
% df = cell(1, Nf);
% for ifile = 1:Nf
%     df{ifile} = load(fullfile(datapath, datafile{ifile}));
% end

% Npixel = 864;
% Nslice = 64;
% Nview = 1152;
% 
% m = 32;
% cut = [0.5e4, 0.5e4, 1.0e4, 1.0e4];
% 
% % bhcorr
% bhcorr = loaddata('E:\matlab\CTsimulation\cali\experiment\testdata\beamharden_Body_120_v1.0.corr');
% p_bh = reshape(double(bhcorr.main), Npixel*Nslice, []);
% A = cell(1, Nf);
% B = cell(1, Nf);
% for ifile = 1:Nf
%     A{ifile} = reshape(df{ifile}.rawdata, Npixel, []);
%     B{ifile} = reshape(df{ifile}.rawdata_bh, Npixel, []);
%     for iv = 1:Nview*Nslice
%         ix_l = find(A{ifile}(:, iv)>cut(ifile), 1, 'first');
%         ix_r = find(A{ifile}(:, iv)>cut(ifile), 1, 'last');
%         A{ifile}([1:ix_l+m-1 ix_r-m+1:Npixel], iv) = 0;
%         B{ifile}([1:ix_l+m-1 ix_r-m+1:Npixel], iv) = 0;
%     end
%     A{ifile} = reshape(A{ifile}, [], Nview);
%     B{ifile} = reshape(B{ifile}, [], Nview);
% end



% non-linear4
midu = 416;
n1 = 2;
n2 = 2;
p_nl = zeros(Npixel*Nslice, n1);
sy = zeros(Npixel*Nslice, Nf);

K = 1.0;
options = optimoptions('lsqnonlin','Display','off');
% options = [];
% ipx = 400;
w = repmat([2 1 2 1]', 1, Nview);

for islice = 1:Nslice
% islice = 6;
    tic
%     ii = 400;
    t0 = zeros(1, n1);
    t0(n1) = 1.0;
    t1 = zeros(1, n2);
    t1(n2) = 1.0;
    for ii = midu+1:Npixel
        ipx = ii+(islice-1)*Npixel;
        x = double([B{1}(ipx, :); B{2}(ipx, :); B{3}(ipx, :); B{4}(ipx, :);]./(1000*K));
        y = double([A{1}(ipx, :); A{2}(ipx, :); A{3}(ipx, :); A{4}(ipx, :);]./(1000*K));
        s = y>0;
        sy(ipx, :) = sum(s, 2);
        if all(sy(ipx, :)<Nview/4)
            continue;
        end
%         p_ipx = lsqnonlin(@(t) (iterpolyval(t, x(s))-y(s)), t0, [], [], options);
        if sum(sy(ipx, :)<Nview/2)<2
            p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), t0, [], [], options);
            if p_ipx(1)>1.0
                p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), [0 0 1], [], [], options);
            end
            t0 = p_ipx;
            p_nl(ipx, :) = p_ipx;
        else
            % order curse
            p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), t1, [], [], options);
            t1 = p_ipx;
            p_nl(ipx, (1+n1-n2:n1)) = p_ipx;
        end
    end
    
    t0 = zeros(1, n1);
    t0(n1) = 1.0;
    t1 = zeros(1, n2);
    t1(n2) = 1.0;
    for ii = midu:-1:1
        ipx = ii+(islice-1)*Npixel;
        x = double([B{1}(ipx, :); B{2}(ipx, :); B{3}(ipx, :); B{4}(ipx, :);]./(1000*K));
        y = double([A{1}(ipx, :); A{2}(ipx, :); A{3}(ipx, :); A{4}(ipx, :);]./(1000*K));
        s = y>0;
        sy(ipx, :) = sum(s, 2);
        if all(sy(ipx, :)<Nview/4)
            continue;
        end
%         p_ipx = lsqnonlin(@(t) (iterpolyval(t, x(s))-y(s)), t0, [], [], options);
        if ii==200
            1;
        end
        if sum(sy(ipx, :)<Nview/2)<2
            p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), t0, [], [], options);
            if p_ipx(1)>1.0
                p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), [0 0 1], [], [], options);
            end
            t0 = p_ipx;
            p_nl(ipx, :) = p_ipx;
        else
            % order curse
            p_ipx = lsqnonlin(@(t) (iterinvpolyval(t, y(s))-x(s)).*w(s), t1, [], [], options);
            t1 = p_ipx;
            p_nl(ipx, (1+n1-n2:n1)) = p_ipx;
        end
    end
    toc;
end

% p_nl(isnan(p_nl)) = 0;
p_nl(p_nl(:,n1)==0, n1) = 1;
