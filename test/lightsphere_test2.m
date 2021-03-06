% sample 0 
Np = 10000;
hset1 = haltonset(6, 'Skip', 0);
hskip1 = 0;
hseq1 = net(hset1, Np);
% hseq1 = [ 0 0 0 1e-9 0 0];

V0 = pi*4/3;
r1 = samplesetinobject(hseq1(:,1:3), 'sphere');
r2 = samplesetinobject(hseq1(:,4:6), 'sphere');
r0 = [2.0 0.0 0];

% #1
J1 = V0^2/Np;
R1 = numexperiment2(r1, r2, r0, J1);

% #2
r1_polar = xyz2polar(r1);
r2_polar = xyz2polar(r2);
[v1_polar, sinp12, J2_p] = S3renormalize(r2_polar, r1_polar);
v1 = polar2xyz(v1_polar);
J2 = sinp12.*J2_p.*(V0^2/Np);
R2 = numexperiment2(v1, r2, r0, J2);

% #3
R3 = numexperiment2_rn(v1, r2, r0, sinp12, J2_p.*(V0^2/Np));

% #3 seq
R3seq = numexperiment2_seq(v1, r2, r0, sinp12, J2_p.*V0^2);

% fun1
function R = numexperiment2(r1, r2, r0, J)

a = sum(r0.^2)*(1/4/pi^2);
r12 = r2 - r1;
r20 = r0 - r2;
R = sum(sum(normr(r20).*normr(r0), 2)./sum(r12.^2,2)./sum(r20.^2,2).*J)*a;

end


% fun2
function R = numexperiment2_rn(r1, r2, r0, sinp, Jp)

a = sum(r0.^2)*(1/4/pi^2);
tol_p = 1e-9;

r12 = r2 - r1;
r20 = r0 - r2;

Rs = sum(normr(r20).*normr(r0), 2)./sum(r12.^2,2)./sum(r20.^2,2).*sinp.*Jp;

s_zero = sinp<tol_p;
if any(s_zero)
    Rs(s_zero) = sum(normr(r20(s_zero,:)).*normr(r0),2).*(8./(1+sum(r2(s_zero,:).^2,2)).^2)./sum(r20(s_zero,:).^2,2).*Jp(s_zero);
end

R = sum(Rs)*a;

end

% fun2seq
function Rs = numexperiment2_seq(r1, r2, r0, sinp, Jp)

a = sum(r0.^2)*(1/4/pi^2);
tol_p = 1e-6;

r12 = r2 - r1;
r20 = r0 - r2;

Rs = sum(normr(r20).*normr(r0), 2)./sum(r12.^2,2)./sum(r20.^2,2).*sinp.*Jp;

s_zero = sinp<tol_p;
if any(s_zero)
    Rs(s_zero) = sum(normr(r20(s_zero,:)).*normr(r0),2).*(8./(1+sum(r2(s_zero,:).^2,2)).^2)./sum(r20(s_zero,:).^2,2).*Jp(s_zero);
end

Rs = Rs.*a;

end