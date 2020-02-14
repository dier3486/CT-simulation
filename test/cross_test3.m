% run in reconnode_crosstalkcali.m
% try3

m_mod = 16;

% air rate
% rawair = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_large_v1.0.raw');
% rawempty = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_empty_v1.0.raw');
airrate = mean([rawair(:).Raw_Data], 2)./mean([rawempty(:).Raw_Data], 2);
airrate = reshape(airrate, Npixel, Nslice);
% for ii = 1:Nslice
%     airrate(:,ii) = smooth(airrate(:,ii), 0.05, 'rloess');
% end
airrate_mg = squeeze(mean(reshape(airrate, Npixel, Nmerge, Nslice_mg), 2));
% airrate_mg = airrate_mg(:);

% diff
crsb1 = [zeros(1, Nslice_mg); 1-airrate_mg(1:end-1,:)./airrate_mg(2:end,:)];
crsb2 = [airrate_mg(2:end,:)./airrate_mg(1:end-1,:)-1; zeros(1, Nslice_mg)];
crsbase = reshape(mean([crsb1(:) crsb2(:)],2), Npixel, Nslice_mg);
% smooth
crsbase(1:m_mod, :) = 0;
crsbase(end-m_mod+1:end, :) = 0;
for ii = 1:Nslice_mg
    crsbase(:,ii) = smooth(crsbase(:,ii), 0.07, 'rloess');
end


alpha = -5.0;

% ini the buffer to use
Nvbk = Nview*Nbk/2;
x = zeros(3, Nvbk);
y = zeros(1, Nvbk);

lambda = 0.02;

start_mod = 18;
end_mod = 19;
Np1 = m_mod*(start_mod-1)+1;
Np2 = Npixel-m_mod*(end_mod-1);
options = optimoptions('lsqnonlin','Display','off');
% options = [];
p_crs = zeros(Npixel, Nslice_mg);
pd = zeros(Npixel, Nslice_mg);
mp = zeros(m_mod, Nslice_mg);
mpd = zeros(m_mod, Nslice_mg);
for islice = 1:Nslice_mg
% for islice = 8
    slice_index = (1:Nmerge) + (islice-1)*Nmerge;
    index_isl = (1:Npixel*Nmerge) + (islice-1)*Npixel*Nmerge;
    
    shift_sl = (islice-1)*Npixel;
    for ipx = Np1:Np2
    % for ipx = 417
        % pixel index
        ipixel = ipx+(islice-1)*Npixel;
        % index range
        Srange = zeros(Npixel, Nview, Nbk/2);
        for ibk = 1:Nbk/2
            for iview = 1:Nview
                Srange(index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk), iview, ibk) = 1;
            end
        end
        % set the data to fit
        for ibk = 1:Nbk/2
            index_xy = (1:Nview) + (ibk-1)*Nview;
            % x is the original data in fitting
            ix = ibk*2-1;   % 1 3 5 7
            x(:, index_xy) = double(dataflow.(datafields{ix})(ipixel-1:ipixel+1, :));
            % y is the target data
            iy = ibk*2;     % 2 4 6 8
            y(index_xy) = double(dataflow.(datafields{iy})(ipixel, :));
            % cut y by index range
            y(index_xy) = y(index_xy).*Srange(ipx, :, ibk);
        end
        % s
        s = y>0;
        x = x(:, s);
        y = y(s);
        p_crs(ipx, islice) = lsqnonlin(@(t) crossfit1(t, x, y, lambda, crsbase(ipx,islice).*alpha), 0, [], [], options);
    end
    p_crs(Np1:Np2, islice) = p_crs(Np1:Np2, islice)-mean(p_crs(Np1:Np2, islice));
    mp(:, islice) = mean(reshape(p_crs(Np1:Np2, islice), m_mod, []), 2);
    
    pd(:, islice) = cumsum(p_crs(:, islice)).*2;
    
end

p0_crs = p_crs;

