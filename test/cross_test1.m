% run in reconnode_crosstalkcali.m

% alpha = -5.0;

m_mod = 16;
% start_mod = 18;
% end_mod = 19;
% Np1 = m_mod*(start_mod-1)+1;
% Np2 = Npixel-m_mod*(end_mod-1);

islice = 2;

index_isl = (1:Npixel) + (islice-1)*Npixel;

Ax = double([dataflow.rawdata_bk1(index_isl, :) dataflow.rawdata_bk3(index_isl, :)]);
Ay = double([dataflow.rawdata_bk2(index_isl, :) dataflow.rawdata_bk4(index_isl, :)]);

Srange = false(Npixel, Nview, Nbk/2);
for ibk = 1:Nbk/2
    for iview = 1:Nview
        Srange(index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk), iview, ibk) = true;
    end
end
Srange = reshape(Srange, Npixel, Nview*Nbk/2);

lambda = 0;
p0 = -ones(Npixel-1, 1).*0.02;
p1 = lsqnonlin(@(t) crossfit2(t, Ax, Ay, Srange, lambda), p0);

