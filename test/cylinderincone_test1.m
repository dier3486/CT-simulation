% An (almost-Z) cylinder cut by collimator blades

Zmod = 5;

SID = 700;
DID = 450;
SDT = 100;

% focal and coll
fpos0 = [0 -SID -SDT.*Zmod];
fangle = 1.0;
Arot = [cos(fangle) sin(fangle) 0; -sin(fangle) cos(fangle) 0; 0 0 1];

fpos = fpos0*Arot;
cone1 = struct();
cone1.type = 'cone';
cone1.O = fpos;
cone1.V = [DID 0    0; 0   DID  0; [0 0 -40.*Zmod]-fpos];
cone1.anglerange = [pi/6 pi*5/6]+fangle;
cone2 = struct();
cone2.type = 'cone';
cone2.O = fpos;
cone2.V = [DID 0    0; 0   DID  0; [0 0 40.*Zmod]-fpos];
cone2.anglerange = [pi/6 pi*5/6]+fangle;

% det
det = struct();
det.type = 'Tube';
det.O = [0 0 0];
det.V = [DID 0   0;
         0   DID 0;
         0   0   -40.*Zmod]; 

% cylinder   
obj1 = struct();
obj1.type = 'Cylinder';
obj1.O = [0 30 -40.*Zmod];
obj1.V = [160  0  0;
      0 160  0;
      0  0  30.*Zmod];
t1 = (rand(1)-0.5).*0.3;
t2 = (rand(1)-0.5).*0.3;
% t2 = 0;
Arot = [cos(t1) 0 sin(t1); 0 1 0; -sin(t1) 0 cos(t1)]*[1 0 0; 0 cos(t2) sin(t2); 0 -sin(t2) cos(t2)];
obj1.V = obj1.V*Arot;

% inv
obj1.invV = inv(obj1.V);

% sample 0 
Np = 5e3;
hset1 = haltonset(3);
hskip1 = 0;
hseq1 = net(hset1, Np);
u0 = samplesetinobject(hseq1, 'spcylinder');
v0 = u0*obj1.V+obj1.O;

% to cone1
zsign = sign(obj1.V(3,3));
t12 = zeros(Np, 2);
% invc1 = inv(cone1);
u1 = (v0-fpos)/cone1.V;
nz1 = [0 0 1]*obj1.V/cone1.V;
a = nz1(1)^2+nz1(2)^2-nz1(3)^2;
b = u1*[nz1(1); nz1(2); -nz1(3)];
c = u1(:,1).^2 + u1(:,2).^2 - u1(:,3).^2;
t12(:, (1-zsign)/2+1) = (-sqrt(b.^2-a.*c).*zsign-b)./a;

u2 = (v0-fpos)/cone2.V;
nz2 = [0 0 1]*obj1.V/cone2.V;
a = nz2(1)^2+nz2(2)^2-nz2(3)^2;
b = u2*[nz2(1); nz2(2); -nz2(3)];
c = u2(:,1).^2 + u2(:,2).^2 - u2(:,3).^2;
t12(:, (1+zsign)/2+1) = (-sqrt(b.^2-a.*c).*zsign-b)./a;

u2 = u0;
t12 = t12+u2(:,3);
t12(t12>1) = 1;
t12(t12<-1) = -1;
u2(:, 3) = ((t12(:,2)-t12(:,1)).*u2(:, 3) + sum(t12,2))./2;
v2 = u2*obj1.V+obj1.O;
w2 = (t12(:,2)-t12(:,1))./2;


% plot
objtoplot = {'det', 'obj1', 'cone1', 'cone2'};
Nobj = length(objtoplot);
f1 = figure; hold on;
map = [];
for ii = 1:Nobj
    [h, map_ii] = objectvisual(eval(objtoplot{ii}), f1);
    h.CData = (h.CData+1).*0.9999 + (ii-1)*2;
    map = [map; map_ii];
end
colormap(map);
caxis([0,2*Nobj]);

plot3(fpos(1), fpos(2), fpos(3), '*');

plot3(v2(:,1),v2(:,2),v2(:,3),'.','MarkerSize', 4);

% axis equal