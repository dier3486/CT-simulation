% sample 0 
Np = 5e3;
hset1 = haltonset(3);
hskip1 = 0;
hseq1 = net(hset1, Np);

V0 = pi*4/3;
u0 = samplesetinobject(hseq1, 'sphere');
r0 = [0.1 0.0 0];

% #1
J1 = V0/Np;
R1 = numexperiment1(u0, r0, J1);

% #2
u0_polar = xyz2polar(u0);
r0_polar = xyz2polar(r0);
[v2_polar, sinp12, J2_p] = S3renormalize(r0_polar, u0_polar);
v2 = polar2xyz(v2_polar);
J2 = sinp12.*J2_p.*(V0/Np);
R2 = numexperiment1(v2, r0, J2);

% #3
R3 = numexperiment1_rn(v2, r0, sinp12, J2_p.*(V0/Np));

% fun1
function R = numexperiment1(r, r0, J)

a = sum(r0.^2)*(3/4/pi);
r1 = r0-r;
R = sum(sum(normr(r1).*normr(r0), 2)./sum(r1.^2,2).*J)*a;

end


% fun2
function R = numexperiment1_rn(r, r0, sinp, Jp)

a = sum(r0.^2)*(3/4/pi);
tol_p = 1e-9;

r1 = r0-r;
Rs = sum(normr(r1).*normr(r0), 2)./sum(r1.^2,2).*sinp.*Jp;

s_zero = sinp<tol_p;
if any(s_zero)
    Rs(s_zero) = sum(normr(r1(s_zero,:)).*normr(r0), 2).*(8/(1+sum(r0.^2))^2).*Jp(s_zero);
end

R = sum(Rs)*a;

end