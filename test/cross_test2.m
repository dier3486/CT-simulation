% run in reconnode_crosstalkcali.m
% try2

% alpha = -5.0;

m_mod = 16;
% start_mod = 18;
% end_mod = 19;
% Np1 = m_mod*(start_mod-1)+1;
% Np2 = Npixel-m_mod*(end_mod-1);

Nmerge = 8;
Nslice_mg = Nslice/Nmerge;

p_crs = zeros(Npixel, Nslice_mg);

for islice = 1:Nslice_mg
    fprintf('.');
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
fprintf('\n');

