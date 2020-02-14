function [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)
% crosstalk calibration
% [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;

% parameters to use
if ~isempty(status)
    caliprm = prmflow.pipe.(status.nodename);
else
    % for debug
    caliprm = struct();     
end

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end

% debug test
Nmerge = 4;
Nslice_mg = Nslice/Nmerge;
% Npixelpermod = 16;

% Nppm = Npixelpermod;

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

% I know the air corr is in prmflow.corrtable
Aircorr = prmflow.corrtable.Air;
Aircorr.main = reshape(Aircorr.main, Nps, []);
airrate = mean(Aircorr.main, 2);
% I know the beamharden corr is in prmflow.corrtable
BHcorr = prmflow.corrtable.Beamharden;
BHcorr.main = reshape(BHcorr.main, Nps, []);
% I know the nonlinear corr is in dataflow
NLcorr = dataflow.nonlinearcorr;
NLcorr.main = reshape(NLcorr.main, Nps, []);
% I know
HCscale = 1000;

% % inverse the ideal data
% for ibk = 2:2:Nbk
%     % inverse Housefield
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
%     for ipx = 1:Nps
%         % inverse nonlinear
%         dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
%         % inverse beamharden
%         dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
%     end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
% end
% % inverse the original data
% for ibk = 1:2:Nbk
%     % inverse Housefield
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
%     for ipx = 1:Nps
%         % inverse beamharden
%         dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
%     end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
% end


% inverse the ideal data
for ibk = 2:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % inverse beamharden
        dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
    end
    % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end
% inverse the original data
for ibk = 1:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % apply the non-linear corr
        dataflow.(datafields{ibk})(ipx, :) = iterpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        % inverse beamharden
        dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
    end
    % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end


% tmp hard codes (from cross_test2)

p_crs = zeros(Npixel, Nslice_mg);

for islice = 1:Nslice_mg
    fprintf(repmat('.', 1, Nmerge));
    slice_index = (1:Nmerge) + (islice-1)*Nmerge;
    
    index_isl = (1:Npixel*Nmerge) + (islice-1)*Npixel*Nmerge;

    x = double([dataflow.rawdata_bk1(index_isl, :) dataflow.rawdata_bk3(index_isl, :)]);
    y = double([dataflow.rawdata_bk2(index_isl, :) dataflow.rawdata_bk4(index_isl, :)]);
    
    % slice merge
    x = reshape(mean(reshape(x, Npixel, Nmerge, []), 2), Npixel, []);
    y = reshape(mean(reshape(y, Npixel, Nmerge, []), 2), Npixel, []);

    % to intensity
    x = 2.^(-x);
    y = 2.^(-y);

    Nvbk = Nview*Nbk/2;
    Srange = zeros(Npixel*Nmerge, Nview, Nbk/2);
    for ibk = 1:Nbk/2
        for iview = 1:Nview
            for isl = 1:Nmerge
                t1 = index_range(1, slice_index(isl), iview, ibk) + Npixel*(isl-1);
                t2 = index_range(2, slice_index(isl), iview, ibk) + Npixel*(isl-1);
                Srange(t1:t2, iview, ibk) = 1;
            end
        end
    end
    Srange = sum(reshape(Srange, Npixel, Nmerge, Nview*Nbk/2), 2) == Nmerge;
    Srange = reshape(Srange, Npixel, Nview*Nbk/2);
    Sedge = sum(Srange, 2)<Nview*Nbk/4;
    Srange(Sedge, :) = 0;

    x = x.*Srange;
    y = y.*Srange;

    weight = sum(Srange, 2)./Nvbk;

    b_d = -[zeros(1, Nvbk); diff(y, 1)];
    b_u = -b_d;
    b_d(1, :) = b_d(2, :);
    b_u(1, :) = 0;
    index_d = (0:Nvbk-1).*Npixel;
    index_u = index_d - 1;
    B = spdiags([b_d b_u], [index_d index_u], Npixel, Npixel*Nvbk);

    lamda = 1e-3;
    BB = B*B';
    BB_lamda = BB + speye(Npixel).*lamda;

    v_yx = (y-x).*Srange;
    Niter_A = 2;
    A = sparse(Npixel, Npixel);
    p_A = zeros(Npixel, Niter_A);
    alpha = 1.0;

    for iter_A = 1:Niter_A
        v1 = v_yx + A*v_yx.*alpha;
        v2 = B*v1(:);
        p = BB_lamda\v2;
        Niter_p = 3;
        for iter_p = 1:Niter_p
            p = p + BB_lamda \ (v2 - BB*p);
        end
        p_A(:, iter_A) = p;
        plow_A = p;
        pup_A = [p(2:end); 0];
        d_A = -plow_A-pup_A;
        A = spdiags([plow_A d_A pup_A], [-1 0 1], Npixel, Npixel);
    end
    
    p_crs(:, islice) = p_A(:, end);
end

p_crs = repelem(p_crs, 1, Nmerge);

% paramters for corr
crosstalkcorr = caliprmforcorr(prmflow, corrversion);
% copy results to corr
crosstalkcorr.Nslice = Nslice;
crosstalkcorr.order = 1;
crosstalkcorr.mainsize = Nps;
crosstalkcorr.main = p_crs;

% to return
dataflow.crosstalkcorr = crosstalkcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end