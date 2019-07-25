% test script 

addpath(genpath('../'));

focal = [0.4, -250, 0];

h = 1;
Xran = [-300, 300];
Zran = [-150, 150];
Xgrid = Xran(1)+h/2 : h : Xran(2)-h/2;
Zgrid = Zran(1)+h/2 : h : Zran(2)-h/2;
Nx = length(Xgrid); Nz = length(Zgrid);
[X, Z] = ndgrid(Xgrid, Zgrid);
Nd = size(X(:),1);
detector = [X(:) ones(Nd,1).*150 Z(:)];

object1.type = 'sphere';

object1.O = [0, 40, 0];
object1.Vector = eye(3).*80;
object1.invV = inv(object1.Vector);
object1.ref_mu = 0.01;
object1.ref_HC = 1000;
object1.material = 'water';

Nv = 1200;
views = linspace(0, pi*2, Nv+1);
views = views(1:end-1);
[D1, L1] = intersection(focal, detector, object1, 'views', views);
