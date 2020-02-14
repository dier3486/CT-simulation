% run in reconnode_crosstalkcali.m

% alpha = -5.0;

m_mod = 16;
% start_mod = 18;
% end_mod = 19;
% Np1 = m_mod*(start_mod-1)+1;
% Np2 = Npixel-m_mod*(end_mod-1);

islice = 2;

index_isl = (1:Npixel) + (islice-1)*Npixel;

x = double([dataflow.rawdata_bk1(index_isl, :) dataflow.rawdata_bk3(index_isl, :)]);
y = double([dataflow.rawdata_bk2(index_isl, :) dataflow.rawdata_bk4(index_isl, :)]);

% to intensity
x = 2.^(-x);
y = 2.^(-y);

Nvbk = Nview*Nbk/2;
Srange = zeros(Npixel, Nview, Nbk/2);
for ibk = 1:Nbk/2
    for iview = 1:Nview
        Srange(index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk), iview, ibk) = 1;
    end
end
Srange = reshape(Srange, Npixel, Nview*Nbk/2);
% Sedge = sum(Srange, 2)<Nview*Nbk/6;
% Srange(Sedge, :) = 0;

weight = sum(Srange, 2)./Nvbk;


% lambda = sqrt(Nview*Nbk) * 1.0;
lambda = 1;
p0 = -ones(Npixel, 1).*0.0;

% p1 = lsqnonlin(@(t) crossfit2(t, Ax, Ay, Srange, lambda), p0);

% try to fit in intensity domain, inverse style

Niter = 10;
p = zeros(Npixel, Niter+1);
p(:, 1) = p0;
alpha = 0.1;

eye_lambda = spdiags(ones(Npixel, 1), 0, Npixel, Npixel);
for ii = 1:Niter
%     p_u = [p(2:end, ii); 0];
%     d = 1 - p(:, ii) - p_u;
%     A = spdiags([p_u d p], [-1 0 1], Npixel, Npixel)';

    p_u = [-p(2:end, ii); 0];
    d = ones(Npixel, 1);
    A = spdiags([p_u d p], [-1 0 1], Npixel, Npixel)';
    
    r_x = A*x - y;
    r_x = r_x.*Srange;
    sum(r_x(:).^2)
    
%     dAy = diff(A\y, 1);
%     x_d = [zeros(1, Nvbk); -dAy];
%     x_u = [zeros(1, Nvbk); dAy];
%     x_d = x_d.*Srange;
%     x_u = x_u.*Srange;
%     index_d = (0:Nvbk-1).*Npixel;
%     index_u = index_d - 1;
%     B = spdiags([x_d x_u], [index_d index_u], Npixel, Npixel*Nvbk)';
    
    Ainv_y = A\y;
    b_d = [zeros(1, Nvbk); (Ainv_y(3:end, :) - Ainv_y(1:end-2, :)); zeros(1, Nvbk)];
    b_d = b_d.*Srange;
    index_d = (0:Nvbk-1).*Npixel;
    B = spdiags(b_d, index_d, Npixel, Npixel*Nvbk)';
    
    BB = B'*B+eye_lambda;
    v = B'*r_x(:) - p(:,ii).*lambda;
    dp = BB\v;
    dp = dp;
    
    p(:, ii+1) = p(:, ii) + dp.*alpha;
end






