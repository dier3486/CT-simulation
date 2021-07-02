% sample set
Nmgc = 61;
Nset = 10;
Nstep = 500;
Np = Nmgc*Nset*Nstep;
hskip0 = 0;

% toy system
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
% detector
det = struct();
det.type = 'Tube';
det.O = [0 0 0];
det.V = [DID 0   0;
         0   DID 0;
         0   0   -40.*Zmod]; 

% object1
% cylinder   
obj1 = struct();
obj1.type = 'Cylinder';
obj1.O = [0 50 -40.*Zmod];
obj1.vector = [120  0  0;
      0 120  0;
      0  0  -30.*Zmod];
t1 = (rand(1)-0.5).*1.0;
t2 = (rand(1)-0.5).*1.0;
% t2 = 0;
Arot = [cos(t1) 0 sin(t1); 0 1 0; -sin(t1) 0 cos(t1)]*[1 0 0; 0 cos(t2) sin(t2); 0 -sin(t2) cos(t2)];
obj1.vector = obj1.vector*Arot;
obj1.invV = inv(obj1.vector);
obj1.volume = defaultvolume(obj1.vector, obj1.type);

% object2
obj2 = obj1;
% object2 = object1

Nr0 = 20;
r0 = [ones(Nr0, 1).*2  linspace(-5, 5, Nr0)'  zeros(Nr0,1)];

R = zeros(Nstep, Nmgc, Nr0);
tic;
for istep = 1:Nstep
    hskip = hskip0 + (istep-1)*Nmgc*Nset;
    hset1 = haltonset(6, 'Skip', hskip);
    hseq1 = net(hset1, Nmgc*Nset);
    % sequence on S3
    phi1 = samplesetinobject(hseq1(:,1:3), 'sphere', true);
%     phi2 = samplesetinobject(hseq1(:,4:6), 'sphere', true);
    % sequence in unit
    u1 = samplesetinobject(hseq1(:,1:3), 'spcylinder', false);
    u2 = samplesetinobject(hseq1(:,4:6), 'spcylinder', false);
    % sequence in objects
    r1 = u1*obj1.vector + obj1.O;
    r2 = u2*obj2.vector + obj2.O;
    
    % phi2=f1(r2);
    phi2 = samplesetinobject((r2 - obj1.O)*obj1.invV, 'cylinder2sphere', true, false);
    
    % renorm
    [phi1_rn, sinp12, Jstr] = S3renormalize(phi2, phi1);
    
    % new r1
    u1_rn = samplesetinobject(phi1_rn, 'sphere2cylinder', false, true);
    r1_rn = u1_rn*obj1.vector + obj1.O;
    
    % J
    J1 = obj1.volume;
    J2 = obj2.volume;
    J = Jstr.*J1.*J2;
    
    % R
    R3seq = numexperiment2_seq(r1_rn, r2, r0, sinp12, J);
    R(istep, :, :) = reshape(sum(reshape(R3seq, Nmgc, Nset, Nr0), 2)./Nset, 1, Nmgc, Nr0);
    if istep>1
        R(istep, :, :) = (R(istep-1, :, :).*(istep-1) + R(istep, :, :))./istep;
    end
    
end
toc;

Rret = squeeze(mean(R, 2));
Rerr = squeeze(std(R,1,2)./Nmgc);



% fun2seq
function Rs = numexperiment2_seq(r1, r2, r0, sinp, Jp)

n1 = size(r1,1);
n0 = size(r0,1);
r12 = r2 - r1;
r20 = repelem(r0, n1, 1) - repmat(r2, n0, 1);

D12 = 1./sum(r12.^2,2).*sinp.*Jp;
D20 = sum(normr(r20).*repelem(normr(r0), n1, 1), 2)./sum(r20.^2,2);

Rs = repmat(D12, n0, 1).*D20;
end