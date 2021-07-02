% sample 0 
Nmgc = 31;
Nset = 10;
Nstep = 1000;
Np = Nmgc*Nset*Nstep;

hskip0 = 0;
V0 = pi*4/3;
r0 = [2.0 0.0 0];

R = zeros(Nstep, Nmgc);
for istep = 1:Nstep
    hskip = hskip0 + (istep-1)*Nmgc*Nset;
    hset1 = haltonset(6, 'Skip', hskip);
    
    hseq1 = net(hset1, Nmgc*Nset);
    r1 = samplesetinobject(hseq1(:,1:3), 'sphere');
    r2 = samplesetinobject(hseq1(:,4:6), 'sphere');
    
    r1_polar = xyz2polar(r1);
    r2_polar = xyz2polar(r2);
    [v1_polar, sinp12, J2_p] = S3renormalize(r2_polar, r1_polar);
    v1 = polar2xyz(v1_polar);
    R3seq = numexperiment2_seq(v1, r2, r0, sinp12, J2_p);
    R(istep, :) = sum(reshape(R3seq, Nmgc, Nset), 2)'.*(V0^2/Nset);
    if istep>1
        R(istep, :) = (R(istep-1, :).*(istep-1) + R(istep, :))./istep;
    end
    
end

Rret = mean(R, 2);
Rerr = std(R,1,2)./Nmgc;



% fun2seq
function Rs = numexperiment2_seq(r1, r2, r0, sinp, Jp)

a = sum(r0.^2)*(1/4/pi^2);
tol_p = 1e-9;

r12 = r2 - r1;
r20 = r0 - r2;

Rs = sum(normr(r20).*normr(r0), 2)./sum(r12.^2,2)./sum(r20.^2,2).*sinp.*Jp;

% s_zero = sinp<tol_p;
% if any(s_zero)
%     Rs(s_zero) = sum(normr(r20(s_zero,:)).*normr(r0),2).*(8./(1+sum(r2(s_zero,:).^2,2)).^2)./sum(r20(s_zero,:).^2,2).*Jp(s_zero);
% end

Rs = Rs.*a;

end