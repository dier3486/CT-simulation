% function [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)
% crosstalk calibration
% [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;

% debug test

Nmerge = 1;
Nslice_mg = Nslice/Nmerge;

% step #2 inverse bh and nl
% I know we should do step1 first, not now

% find out the input data, I know they are rawdata_bk1, rawdata_bk2 ...
datafields = findfields(dataflow, '\<rawdata_bk');
Nbk = length(datafields);
headfields = findfields(dataflow, '\<rawhead_bk');

% in which odd are the original value, even are the ideal value (fitting target)
% reshape
for ibk = 1:Nbk
    dataflow.(datafields{ibk}) = reshape(dataflow.(datafields{ibk}), Nps, Nview);
end
% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

% I know the beamharden corr is in prmflow.corrtable
BHcorr = prmflow.corrtable.Beamharden;
BHcorr.main = reshape(BHcorr.main, Nps, []);
% I know the nonlinear corr is in dataflow
NLcorr = dataflow.nonlinearcorr;
NLcorr.main = reshape(NLcorr.main, Nps, []);
% I know
HCscale = 1000;
% inverse the ideal data
for ibk = 2:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % inverse nonlinear
        dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        % inverse beamharden
        dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
    end
end
% inverse the original data
for ibk = 1:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % inverse beamharden
        dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
    end
end


% step #1 air rate
% I know
rawair = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_large_v1.0.raw');
rawempty = loaddata('E:\data\rawdata\bhtest\rawdata_staticair_120KV200mA_empty_v1.0.raw');
airrate = (mean([rawair(:).Raw_Data], 2)-16384)./(mean([rawempty(:).Raw_Data], 2)-16384);
airrate = reshape(airrate, Npixel, Nslice);
% for ii = 1:Nslice
%     airrate(:,ii) = smooth(airrate(:,ii), 0.05, 'rloess');
% end
crsb1 = [zeros(1, Nslice); 1-airrate(1:end-1,:)./airrate(2:end,:)];
crsb2 = [airrate(2:end,:)./airrate(1:end-1,:)-1; zeros(1, Nslice)];
crsbase = reshape(mean([crsb1(:) crsb2(:)],2), Npixel, Nslice);
for ii = 1:Nslice
    crsbase(:,ii) = smooth(crsbase(:,ii), 0.05, 'rloess');
end


% those codes will be moved to beamharden calibration and the airrate shall be saved in beamharden calibration table

% alpha = -5.0;
% lambda = 0.02;
% m_mod = 16;
% start_mod = 18;
% end_mod = 19;
% Np1 = m_mod*(start_mod-1)+1;
% Np2 = Npixel-m_mod*(end_mod-1);
% 
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


% % status
% status.jobdone = true;
% status.errorcode = 0;
% status.errormsg = [];
% end