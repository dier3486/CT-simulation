% nonlinear calibration test

% load data
datapath = 'E:\data\rawdata\bhtest\';
datafile = {'df_w200c.mat', 'df_w200off100.mat', 'df_w300c.mat', 'df_w300off100.mat'};
Nf = length(datafile);
df = cell(1, Nf);
for ifile = 1:Nf
    df{ifile} = load(fullfile(datapath, datafile{ifile}));
end

Npixel = 864;
Nslice = 64;
Nview = 1152;

m = 16;
cut = 1e3;

% % manual bad channel
% badchannel = [2919 12609];
% badmod28 = 448+(32:47).*864;
% badmod19 = 304+(48:63).*864;
% badchannel = [badchannel badmod28 badmod19];
badchannel = [];

A = cell(1, Nf);
B = cell(1, Nf);
for ifile = 1:Nf
    A{ifile} = reshape(df{ifile}.rawdata, Npixel, []);
    B{ifile} = reshape(df{ifile}.rawdata_bh, Npixel, []);
    for iv = 1:Nview*Nslice
        ix_l = find(A{ifile}(:, iv)>cut, 1, 'first');
        ix_r = find(A{ifile}(:, iv)>cut, 1, 'last');
        A{ifile}([1:ix_l+m-1 ix_r-m+1:Npixel], iv) = 0;
        B{ifile}([1:ix_l+m-1 ix_r-m+1:Npixel], iv) = 0;
    end
    A{ifile} = reshape(A{ifile}, [], Nview);
    B{ifile} = reshape(B{ifile}, [], Nview);
    % fix bad channel (tmp)
    for ii = 1:length(badchannel)
        bad_ii = badchannel(ii);
        B{ifile}(bad_ii,:) = (B{ifile}(bad_ii-1,:)+B{ifile}(bad_ii+1,:))./2;
    end
end


% % non-linear1
% n = 2;
% p_nl = nan(Npixel*Nslice, n);
% t0 = zeros(1, n);
% t0(n) = 1.0;
% K = 1.0;
% options = optimoptions('lsqnonlin','Display','off');
% % ipx = 400;
% sy = zeros(Npixel,1);
% for ipx = 1:Npixel*Nslice
%     x = double([B{2}(ipx, :) B{4}(ipx, :)]./(1000*K));
%     y = double([A{2}(ipx, :) A{4}(ipx, :)]./(1000*K));
%     s = y>0;
%     sy(ipx) = sum(s);
%     if sy(ipx)<Nview/2
%         continue;
%     end
%     p_ipx = lsqnonlin(@(t) (iterpolyval(t, x(s))-y(s)), t0, [], [], options);
% %     t0 = p_ipx;
%     p_nl(ipx, :) = p_ipx;
% end
% 
% p_nl = reshape(fillmissing(reshape(p_nl, Npixel, []), 'nearest'), [], n);
% 
% % % or load saved p
% % pdata = load('E:\data\rawdata\bhtest\p2_n2.mat');
% % p_nl = pdata.p;
% 
% % bhcorr
% bhcorr = loaddata('E:\matlab\CTsimulation\cali\experiment\testdata\beamharden_Body_120_v1.0.corr');
% p_bh = reshape(double(bhcorr.main), Npixel*Nslice, []);
% 
% p_scale = p_bh(:,end).*p_nl(:,end);
% 
% % marged A B
% Nmerge = 8;
% Nslice_mg = Nslice/Nmerge;
% Amg = cell(1, Nf);
% Bmg = cell(1, Nf);
% for ifile = [2, 4]
%     % Amerge
%     Amg{ifile} = A{ifile}./1000./p_scale;
%     Amg{ifile}(A{ifile}==0) = nan;
%     Amg{ifile} = reshape(mean(reshape(Amg{ifile}, Npixel, Nmerge, Nslice_mg*Nview), 2), Npixel*Nslice_mg, Nview);
%     Amg{ifile}(isnan(Amg{ifile})) = 0;
%     % Bmerge
%     Bmg{ifile} = B{ifile}./1000;
%     Bmg{ifile}(B{ifile}==0) = nan;
%     % correct B
%     Bmg{ifile} = iterpolyval(p_nl, Bmg{ifile});
%     Bmg{ifile} = Bmg{ifile}./p_scale;
%     Bmg{ifile} = reshape(mean(reshape(Bmg{ifile}, Npixel, Nmerge, Nslice_mg*Nview), 2), Npixel*Nslice_mg, Nview);
%     Bmg{ifile}(isnan(Bmg{ifile})) = 0;
% end
% 
% 
% % cross talk
% 
% % air rate
% rawair = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_large_v1.0.raw');
% rawempty = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_empty_v1.0.raw');
% airrate = mean([rawair(:).Raw_Data], 2)./mean([rawempty(:).Raw_Data], 2);
% airrate = reshape(airrate, Npixel, Nslice);
% for ii = 1:Nslice
%     airrate(:,ii) = smooth(airrate(:,ii), 0.05, 'rloess');
% end
% 
% airrate_mg = squeeze(mean(reshape(airrate, Npixel, Nmerge, Nslice_mg), 2));
% % airrate_mg = airrate_mg(:);
% crsb1 = [zeros(1, Nslice_mg); 1-airrate_mg(1:end-1,:)./airrate_mg(2:end,:)];
% crsb2 = [airrate_mg(2:end,:)./airrate_mg(1:end-1,:)-1; zeros(1, Nslice_mg)];
% crsbase = reshape(mean([crsb1(:) crsb2(:)],2), Npixel, Nslice_mg);
% alpha = -5.0;
% 
% lambda = 0.02;
% m_mod = 16;
% start_mod = 18;
% end_mod = 19;
% Np1 = m_mod*(start_mod-1)+1;
% Np2 = Npixel-m_mod*(end_mod-1);
% % options = optimoptions('lsqnonlin','Display','off');
% options = [];
% p_crs = zeros(Npixel, Nslice_mg);
% mp = zeros(m_mod, Nslice_mg);
% for islice = 1:Nslice_mg
% % for islice = 8
%     shift_sl = (islice-1)*Npixel;
%     for ipx = Np1:Np2
%     % for ipx = 417
%         x = double([Bmg{2}(shift_sl+ipx-1:shift_sl+ipx+1, :) Bmg{4}(shift_sl+ipx-1:shift_sl+ipx+1, :)]);
%         y = double([Amg{2}(shift_sl+ipx, :) Amg{4}(shift_sl+ipx, :)]);
%         s = all(x>0);
%         x = x(:, s);
%         y = y(s);
%         p_crs(ipx, islice) = lsqnonlin(@(t) crossfit1(t, x, y, lambda, crsbase(ipx,islice).*alpha), 0, [], [], options);
%     end
%     p_crs(Np1:Np2, islice) = p_crs(Np1:Np2, islice)-mean(p_crs(Np1:Np2, islice));
%     mp(:, islice) = mean(reshape(p_crs(Np1:Np2, islice), m_mod, []), 2);
% end
% 
% p0_crs = p_crs;
% 
% p_crs(1:Np1-1,:) = repmat(mp, start_mod-1, 1);
% p_crs(Np2+1:end,:) = repmat(mp, end_mod-1, 1);
% p_crs = p_crs+crsbase.*alpha;
% p_crs = repelem(p_crs, 1, Nmerge);


