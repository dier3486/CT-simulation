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
cone(1).V = [DID 0    0; 0   DID  0; [0 0 -40.*Zmod]-fpos0];
cone(1).anglerange = [pi/6+0.05 pi*5/6-0.05];
cone(2) = cone(1);
cone(2).V = [DID 0    0; 0   DID  0; [0 0 40.*Zmod]-fpos0];

% rotation focal and coll
Arot = [cos(fangle) sin(fangle) 0; -sin(fangle) cos(fangle) 0; 0 0 1];
fpos = fpos0*Arot;
cone = objectrotation(cone, Arot);

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
      0  0  -30.*Zmod];
t1 = (rand(1)-0.5).*3;
t2 = (rand(1)-0.5).*3;
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

% u0 cut
% v0 = u0*obj1.V+obj1.O;
% % Z direction
% zsign = sign(obj1.V(3,3));
% if zsign<0
%     cone = cone([2 1]);
% end
% 
% % tcone
% tcone = zeros(Np, 2);
% for ic = 1:2
%     u1 = (v0-fpos)/cone(ic).V;
%     nz1 = [0 0 1]*obj1.V/cone(ic).V;
%     
%     a = nz1(1)^2+nz1(2)^2-nz1(3)^2;
%     b = u1*[nz1(1); nz1(2); -nz1(3)];
%     c = u1(:,1).^2 + u1(:,2).^2 - u1(:,3).^2;
%     d = real(sqrt(b.^2-a.*c));
%     t1 = (-d-b)./a;
%     t2 = (d-b)./a;
%     k1 = u1 + t1*nz1;
%     k2 = u1 + t2*nz1;
%     s1 = k1(:, 3)>=0 & k1(:, 2)>=0 & d>0;
%     s2 = k2(:, 3)>=0 & k2(:, 2)>=0 & d>0;
%     if zsign>=0
%         t1(~s1) = inf;
%         t2(~s2) = inf;
%         tcone(:, ic) = min(t1, t2);
%     else
%         t1(~s1) = -inf;
%         t2(~s2) = -inf;
%         tcone(:,ic) = max(t1, t2);
%     end
% end
% tcone = tcone + u0(:, 3);
% tcone(tcone>1) = 1;
% tcone(tcone<-1) = -1;
% 
% % u2 v2 w2
% u2 = u0;
% u2(:, 3) = ((tcone(:,2)-tcone(:,1)).*u2(:, 3) + sum(tcone,2))./2; 
% w2 = (tcone(:,2)-tcone(:,1))./2;

[u2, w2] = cylinderconcut1(u0, obj1, cone);

v2 = u2*obj1.V+obj1.O;

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
plot3(v2(:,1),v2(:,2),v2(:,3),'.','MarkerSize', 4);
