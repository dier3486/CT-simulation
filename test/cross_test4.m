% run in reconnode_crosstalkcali.m
% try4

% I know Npixelpermod = 16;
Nmod = Npixel/Npixelpermod;

% air rate
airrate = double(prmflow.corrtable.Beamharden.airrate);
airrate = reshape(2.^(-airrate), Npixel, Nslice);
% useless

% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

% loop slice
Pcrs = zeros(Npixel*Nslice, 2);
islice = 2;
for islice = 1:Nslice
    % Seff
    Seff = false(Npixel, Nview*Nbk/2);
    for ibk = 1:Nbk/2
        for iview=1:Nview
            viewindex = iview+(ibk-1)*Nview;
            index_ii = index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk);
            Seff(index_ii, viewindex) = true;
        end
    end

%     Ap = spdiags([p(:,2) -sum(p,2) p(:,1)], [-1 0 1], Npixel,Npixel)';
    
    p = nan(Npixel, 2);
    lamda = 0.0;
    options = optimoptions('lsqnonlin','Display','off');
    % loop pixel
    for ipixel = 2:Npixel-1
        x = zeros(3, Nview*Nbk/2);
        y = zeros(3, Nview*Nbk/2);
        pixelindex = ipixel-1:ipixel+1;
        for ibk = 1:Nbk/2
            ix = ibk*2-1;
            iy = ibk*2;
            viewindex = (1:Nview) + (ibk-1)*Nview;
            x(:, viewindex) = double(dataflow.(datafields{ix})(pixelindex + Npixel*(islice-1), :));
            y(:, viewindex) = double(dataflow.(datafields{iy})(pixelindex + Npixel*(islice-1), :));
        end
    %     rrate = airrate(pixelindex, islice);
        rrate = ones(3,1);
        s = all(Seff(pixelindex, :), 1);
        if any(sum(reshape(s, Nbk/2, Nview),2)<Nview/2)
            continue;
        end
        t0 = [0, 0];
        p(ipixel, :) = lsqnonlin(@(t) crossfit3(t, x, y, s, rrate, lamda), t0, [], [], options);
    end
    
%     p = p - mean(sum(p,2), 'omitnan')/2;
    p = reshape(p, Npixelpermod, Nmod, 2);
    s_unv = any(any(isnan(p),1),3);
    p(:, s_unv, :) = nan;
    pmod = mean(p, 2, 'omitnan');
    p(:, s_unv, :) = repmat(pmod, 1, sum(s_unv));
    p = reshape(p, Npixel, 2);
    p(1,1) = 0; p(end,2) = 0;
    
    pindex = (1:Npixel) + (islice-1)*Npixel;
    Pcrs(pindex, :) = reshape(p, Npixel, 2);
end

% slice merge
Nmerge = 8;
Pcrs_mg = reshape(Pcrs, Npixel, Nmerge, []);
Pcrs_mg = reshape(repmat(mean(Pcrs_mg, 2), 1, Nmerge), [], 2);

% save
tmp = struct();
tmp.pcrs = Pcrs_mg;
save('E:\matlab\CT\SINO\PG\calibration\Pcrs_try4.mat', '-struct', 'tmp');
