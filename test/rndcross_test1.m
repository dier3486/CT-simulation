Npixel = 912;
Nslice = 24;
Npperm = 16;
Nmod = Npixel/Npperm;
Nps = Npixel*Nslice;
Nslamg = 4;

c1 = 0.1;
c2 = 0.05;
p1 = reshape(ones(Nps, 2).*c1, Npixel, Nslice, 2);
p1(Npperm:Npperm:end, :, 1) = 0;
p1(1:Npperm:end, :, 2) = 0;
p2 = randn(Npixel, Nslice, 2).*c2;
p2(1,:,:) = 0;
p2(end,:,:) = 0;
p = reshape(p1+p2, Nps, 2);

Ac = spdiags([-p(:,1) ones(Nps, 1) -p(:,2)], [-1 0 1], Nps, Nps)';
Ad = spdiags(1./(1-sum(p,2)), 0, Nps, Nps);

A = Ac*Ad;
crs = struct();
crs.crossmatrix = A;