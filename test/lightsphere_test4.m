% sample 0 
Nmgc = 61;
Nset = 10;
Nstep = 500;
Np = Nmgc*Nset*Nstep;

hskip0 = 0;
V0 = pi*4/3;
Nr0 = 20;
r0 = [ones(Nr0, 1).*2  linspace(-5, 5, Nr0)'  zeros(Nr0,1)];

R = zeros(Nstep, Nmgc, Nr0);
tic;
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
    R(istep, :, :) = reshape(sum(reshape(R3seq, Nmgc, Nset, Nr0), 2).*(V0^2/Nset), 1, Nmgc, Nr0);
    if istep>1
        R(istep, :, :) = (R(istep-1, :, :).*(istep-1) + R(istep, :, :))./istep;
    end
    
end
toc;

Rret = squeeze(mean(R, 2));
Rerr = squeeze(std(R,1,2)./Nmgc);



% fun2seq
function Rs = numexperiment2_seq(r1, r2, r0, sinp, Jp)

a = sum(r0.^2, 2).*(1/4/pi^2);
tol_p = 1e-9;

n1 = size(r1,1);
n0 = size(r0,1);
r12 = r2 - r1;
r20 = repelem(r0, n1, 1) - repmat(r2, n0, 1);

D12 = 1./sum(r12.^2,2).*sinp.*Jp;
% s_zero = sinp<tol_p;
% if any(s_zero)
%     D12(s_zero) = (8./(1+sum(r2(s_zero,:).^2,2)).^2).*Jp(s_zero);
% end

D20 = sum(normr(r20).*repelem(normr(r0), n1, 1), 2)./sum(r20.^2,2);

Rs = repmat(D12, n0, 1).*D20.*repelem(a, n1, 1);
% Rs = repmat(D12, n0, 1).*D20;


end