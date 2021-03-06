% % sample set
% Nmgc = 61;
% Nset = 10;
% Nstep = 500;
% Np = Nmgc*Nset*Nstep;
% hskip0 = 0;
% 
% % toy system
% Zmod = 5;
% SID = 700;
% DID = 450;
% SDT = 100;
% % focal and coll
% focalposition = [0 -SID -SDT.*Zmod];
% fangle = 1.0;
% clear cone
% cone(2) = struct();
% cone(1).type = 'cone';
% cone(1).O = focalposition;
% cone(1).vector = [DID 0    0; 0   DID  0; [0 0 -40.*Zmod]-focalposition];
% cone(1).anglerange = [pi/6+0.05 pi*5/6-0.05];
% cone(2) = cone(1);
% cone(2).vector = [DID 0    0; 0   DID  0; [0 0 40.*Zmod]-focalposition];
% % rotation focal and coll
% Arot = [cos(fangle) sin(fangle) 0; -sin(fangle) cos(fangle) 0; 0 0 1];
% fpos = focalposition*Arot;
% cone = objectrotation(cone, Arot);
% % detector
% det = struct();
% det.type = 'Tube';
% det.O = [0 0 0];
% det.vector = [DID 0   0;
%          0   DID 0;
%          0   0   40.*Zmod]; 
%     
% % object1
% % cylinder   
% obj1 = struct();
% obj1.type = 'Cylinder';
% obj1.O = [0 50 -40.*Zmod];
% obj1.vector = [120  0  0;
%       0 120  0;
%       0  0  30.*Zmod];
% obj1.Cimage = [];
% t1 = (rand(1)-0.5).*1.0;
% t2 = (rand(1)-0.5).*1.0;
% % t2 = 0;
% Arot = [cos(t1) 0 sin(t1); 0 1 0; -sin(t1) 0 cos(t1)]*[1 0 0; 0 cos(t2) sin(t2); 0 -sin(t2) cos(t2)];
% obj1.vector = obj1.vector*Arot;
% obj1.invV = inv(obj1.vector);
% obj1.volume = defaultvolume(obj1.vector, obj1.type);
% 
% % object2
% obj2 = obj1;
%  
% % plot
% objtoplot = {'det', 'obj1', 'cone(1)', 'cone(2)'};
% Nobj = length(objtoplot);
% f1 = figure; hold on;
% map = [];
% for ii = 1:Nobj
%     [h, map_ii] = objectvisual(eval(objtoplot{ii}), f1);
%     h.CData = (h.CData+1).*0.9999 + (ii-1)*2;
%     map = [map; map_ii];
% end
% colormap(map);
% caxis([0,2*Nobj]);
% % --

[Xr, Yr, Zr] = objectmeshgrid(det, 4, 24);
netsize = size(Xr);
rdet = [Xr(:) Yr(:) Zr(:)];
Ndet = size(rdet, 1);
% -Vn_det
normV_det = normr([Xr(:) Yr(:) zeros(Ndet, 1)]);

Nscatter = 4;
R = cell(Nscatter, 1);
R(:) = {zeros(Nstep, Nmgc, Ndet)};
tol = ones(1, Nscatter).*0;
% K = 0.005;
% Wvol = (obj1.volume.*K).^(1:Nscatter);

Lref = 1000;
mu_ref1 = 0.02;
mu_ref2 = 0.01;

Wvol = (obj1.volume.*mu_ref2).^(1:Nscatter);
% scatter 1

% scatter 2-n

% ini
u1 = cell(1, Nscatter);
u1_inv = cell(1, Nscatter-1);


rc = cell(1, Nscatter);
sinp1 = cell(1, Nscatter-2);
Jstr1 = cell(1, Nscatter-2);

u2 = cell(1, Nscatter);

sinp2 = cell(1, Nscatter);
Jstr2 = cell(1, Nscatter);
w_cone = cell(1, Nscatter);

tol1 = 0.01;
tol2 = 0.2e-3;
% to11 = 0;
% tol2 = 0;
err1 = zeros(Nstep, Nscatter);
err2 = zeros(Nstep, Nscatter);
Sact = true(1, Nscatter);

tic;
for istep = 1:Nstep
    hskip = hskip0 + (istep-1)*Nmgc*Nset;
    hset1 = haltonset(3*Nscatter, 'Skip', hskip);
    hseq1 = net(hset1, Nmgc*Nset);
    
    % Nactsct
    Nactsct = find(Sact,1,'last');
    % ini
    phi1 = cell(1, Nactsct);
    r1 = cell(1, Nactsct-1);
    r2 = cell(1, Nactsct);
    
    % sequence 1 in unit
    u1 = samplesetinobject(hseq1(:,1:3), 'spcylinder', false);
    % sequence on S3
    for isct = 1:Nactsct
        phi1{isct} = samplesetinobject(hseq1(:,(1:3)+(isct-1)*3), 'sphere', true);
    end
    phi2 = phi1;
    
    % r1
    if any(Sact(2:end))
        r1{1} = u1*obj1.vector + obj1.O;
%         for isct = 1:Nscatter-2
        for isct = 1:Nactsct-2
            % I know u2 = (r1-obj1.O)*obj1.invV = u1, so phi2 = phi1
            1;
            % renorm in cylinder
            [phi1{isct+1}, sinp1{isct}, Jstr1{isct}] = S3renormalize(phi1{isct}, phi1{isct+1});
            % phi to r(xyz)
            r1{isct+1} = samplesetinobject(phi1{isct+1}, 'sphere2cylinder', false, true)*obj1.vector + obj1.O;
        end
    end
    
    % r2{1}
    if Sact(1)
        [u2, w_cone{1}] = cylinderconecut(u1, obj1, cone, true);
        r2{1} = u2*obj1.vector + obj1.O;
    end
    % r2{2-m}
%     for isct = 2:Nscatter
    for isct = find(Sact(2:end))+1
        % cone inverse
        u1_inv = cylinderconecut((r1{isct-1} - obj1.O)*obj1.invV, obj1, cone, true);
        phi1_inv = samplesetinobject(u1_inv, 'cylinder2sphere', true, false);
        % renorm in conecut-cylinder
        [phi2{isct}, sinp2{isct}, Jstr2{isct}] = S3renormalize(phi1_inv, phi2{isct});
        % phi to r(xyz)
        [u2, w_cone{isct}] = cylinderconecut(samplesetinobject(phi2{isct}, 'sphere2cylinder', false, true), obj1, cone);
        r2{isct} = u2*obj1.vector + obj1.O;
    end
    
    % Rs
    if Sact(1)
        [Rs_det2, D_det2] = numexperiment2_rdet(r2{1}, rdet, normV_det, obj1);
    end
    if any(Sact(2:end))
        [Rs_det1, D_det1] = numexperiment2_rdet(r1{1}, rdet, normV_det, obj1);
    end
    
    [Rs_rrn, D_rrn] = numexperiment2_rrn(r1, sinp1, Jstr1, obj1);
    [Rs_r12, D_r12] = numexperiment2_r12(r1, r2, sinp2, Jstr2, obj1, Sact);
    [Rs_rof, D_rof] = numexperiment2_rof(r2, focalposition, w_cone, obj1, Sact);
    
    % 1st scatter
    if Sact(1)
        R{1}(istep, :, :) = reshape(sum(reshape(Rs_det2.*repmat(Rs_rof{1}, Ndet, 1), Nmgc, Nset, Ndet), 2).*(Wvol(1)/Nset), ...
            1, Nmgc, Ndet);
        R{1}(istep, :, :) = R{1}(istep, :, :).*reshape(sum(reshape(exp(-mu_ref1.*(D_det2+D_rof{1})), Nmgc, Nset, Ndet), 2), ...
            1, Nmgc, Ndet);
    end
    % 2-m scatter
%     for isct = 2:Nscatter
    for isct = find(Sact(2:end))+1
        Rtmp = Rs_rrn{isct}.*Rs_r12{isct}.*Rs_rof{isct};
%         Rtmp = Rs_r12{isct}.*Rs_rof{isct};
        R{isct}(istep, :, :) = reshape(sum(reshape(Rs_det1.*repmat(Rtmp, Ndet, 1), Nmgc, Nset, Ndet), 2).*(Wvol(isct)/Nset), ...
                               1, Nmgc, Ndet);
        Dtmp = D_det1+D_rrn{isct}+D_r12{isct}+D_rof{isct};
        R{isct}(istep, :, :) = R{isct}(istep, :, :).*reshape(sum(reshape(exp(-mu_ref1.*Dtmp), Nmgc, Nset, Ndet), 2), ...
                               1, Nmgc, Ndet);
    end
    
    % cum mean
    if istep>1
        for isct = 1:Nscatter
            if Sact(isct)
                R{isct}(istep, :, :) = (R{isct}(istep-1, :, :).*(istep-1) + R{isct}(istep, :, :))./istep;
            else
                R{isct}(istep, :, :) = R{isct}(istep-1, :, :);
            end
        end
    end
    
    % error estimate
    Rall = 0;
    for isct = 1:Nscatter
        Rall = Rall + squeeze(mean(R{isct}(istep, :, :), 2));   
    end
    for isct = 1:Nscatter
        if Sact(isct)
            err_isct = squeeze(std(R{isct}(istep, :, :),1,2)./Nmgc);
            err1(istep, isct) = max(err_isct./squeeze(mean(R{isct}(istep, :, :), 2)));
            err2(istep, isct) = max(err_isct./Rall);
        else
            err1(istep, isct) = err1(istep-1, isct);
            err2(istep, isct) = err2(istep-1, isct);
        end
    end
    Sact = Sact & ((err1(istep, :)>tol1) | (err2(istep, :)>tol2));
    if ~any(Sact)
        break;
    end
end
toc;

% result analysis
Rret = cell(1, Nscatter);
Rerr = cell(1, Nscatter);
Intsct = 0;
for isct = 1:Nscatter
    Rret{isct} = squeeze(mean(R{isct}, 2)).*(Lref^2*pi*4);
    Rerr{isct} = squeeze(std(R{isct},1,2)./Nmgc).*(Lref^2*pi*4);
    
    Intsct = Intsct + Rret{isct};
end
Intsct_end = reshape(Intsct(istep, :), netsize);

Rrlt1 = cell(1, Nscatter);
Rrlt2 = cell(1, Nscatter);
for isct = 1:Nscatter
    Rrlt1{isct} = abs(Rerr{isct}./Rret{isct});
    Rrlt2{isct} = Rerr{isct}./Intsct;
end

% model function(s)
% fun1
function [Rs_det1, D_det1] = numexperiment2_rdet(r1, rdet, normV, object)
% r_d^*n^/r_d^2

n0 = size(rdet, 1);
n1 = size(r1,1);

[D_det1, L1] = intersectionABO(rdet, r1, object, 'net');

r10 = normr(repelem(rdet, n1, 1) - repmat(r1, n0, 1));
Rs_det1 = sum(r10.*repelem(normV, n1, 1), 2)./L1.^2./(pi*4);

end
% fun2
function [Rs_rrn, D_rrn] = numexperiment2_rrn(r1, sinp, J, object)
% prod(1/r^2)

% I know r1 is a cell
Nsct = length(r1)+1;

Rs_rrn = cell(1, Nsct);
D_rrn = cell(1, Nsct);
Rs_rrn{1} = 1;
Rs_rrn{2} = 1;
D_rrn{1} = 0;
D_rrn{2} = 0;

for isct = 3:Nsct
%     r12 = r1{isct-2} - r1{isct-1};
%     Rs_rrn{isct} = Rs_rrn{isct-1}./sum(r12.^2, 2).*sinp{isct-2}.*J{isct-2}./(pi*4);
    [Dtmp, L] = intersectionABO(r1{isct-2}, r1{isct-1}, object, 'lines');
    D_rrn{isct} = D_rrn{isct-1} + Dtmp;
    Rs_rrn{isct} = Rs_rrn{isct-1}./L.^2.*sinp{isct-2}.*J{isct-2}./(pi*4);
    
end

end
% fun3
function [Rs_r12, D_r12] = numexperiment2_r12(r1, r2, sinp, J, object, Sact)
% 1/r12^2, single view

Nsct = length(r2);

Rs_r12 = cell(1, Nsct);
D_r12 = cell(1, Nsct);
Rs_r12{1} = 1;
D_r12{1} = 0;

% for isct = 2:Nsct
for isct = find(Sact(2:end))+1
%     r12 = r1{isct-1} - r2{isct};
%     Rs_r12{isct} = 1./sum(r12.^2, 2).*sinp{isct}.*J{isct}./(pi*4);
    [D_r12{isct}, L] = intersectionABO(r1{isct-1}, r2{isct}, object, 'lines');
    Rs_r12{isct} = 1./L.^2.*sinp{isct}.*J{isct}./(pi*4);
end

end
% fun4
function [Rs_rof, D_rof] = numexperiment2_rof(r2, focalpos, w, object, Sact)
% 1/r_f^2, single view

Nsct = length(r2);

Rs_rof = cell(1, Nsct);
D_rof = cell(1, Nsct);
% for isct = 1:Nsct
for isct = find(Sact)
%     rf = r2{isct} - focalpos;
%     Rs_rof{isct} = 1./sum(rf.^2, 2).*w{isct}./(pi*4);
    [D_rof{isct}, L] = intersectionABO(focalpos, r2{isct}, object, 'lines');
    Rs_rof{isct} = 1./L.^2.*w{isct}./(pi*4);
end
end