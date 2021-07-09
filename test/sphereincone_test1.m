% An (almost-Z) cylinder cut by collimator blades

Zmod = 5;

SID = 700;
DID = 450;
SDT = 100;

% focal and coll
fpos0 = [0 -SID -SDT.*Zmod];
fangle = 1.0;
clear cone
cone(2) = struct();
cone(1).type = 'cone';
cone(1).O = fpos0;
cone(1).vector = [DID 0    0; 0   DID  0; [0 0 -40.*Zmod]-fpos0];
cone(1).anglerange = [pi/6+0.05 pi*5/6-0.05];
cone(2) = cone(1);
cone(2).vector = [DID 0    0; 0   DID  0; [0 0 40.*Zmod]-fpos0];

% rotation focal and coll
Arot = [cos(fangle) sin(fangle) 0; -sin(fangle) cos(fangle) 0; 0 0 1];
fpos = fpos0*Arot;
cone = objectrotation(cone, Arot);

% det
det = struct();
det.type = 'Tube';
det.O = [0 0 0];
det.vector = [DID 0   0;
         0   DID 0;
         0   0   -40.*Zmod]; 

% cylinder   
obj1 = struct();
obj1.type = 'Sphere';
obj1.O = [0 30 -40.*Zmod];
obj1.vector = [200  0  0;
      0 200  0;
      0  0  -30.*Zmod];
t1 = (rand(1)-0.5);
t2 = (rand(1)-0.5);
% t2 = 0;
Arot = [cos(t1) 0 sin(t1); 0 1 0; -sin(t1) 0 cos(t1)]*[1 0 0; 0 cos(t2) sin(t2); 0 -sin(t2) cos(t2)];
obj1.vector = obj1.vector*Arot;

% inv
obj1.invV = inv(obj1.vector);

% sample 0 
Np = 5e3;
hset1 = haltonset(3);
hskip1 = 0;
hseq1 = net(hset1, Np);
u0_p = samplesetinobject(hseq1, 'sphere', 1);

% TBC

% i don't like it


% % u0 = u0.*2.0;
% [u2, w2] = cylinderconecut(u0, obj1, cone);
% u3 = cylinderconecut(u2, obj1, cone, 1);
% u4 = cylinderconecut(u0, obj1, cone, 1);

u0 = polar2xyz(u0_p);
v0 = u0*obj1.vector+obj1.O;
% v2 = u2*obj1.vector+obj1.O;
% v3 = u3*obj1.vector+obj1.O;
% v4 = u4*obj1.vector+obj1.O;

% plot
objtoplot = {'det', 'obj1', 'cone(1)', 'cone(2)'};
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
plot3(v0(:,1),v0(:,2),v0(:,3),'.');
% plot3(v3(:,1),v3(:,2),v3(:,3),'.','MarkerSize', 3);
