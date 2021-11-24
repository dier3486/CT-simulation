% 
samplekeV = (5:0.5:140)';
KV = 120;

material_Si = materialdefine(loadmaterial('crystalSI'), samplekeV);
material_Al = materialdefine(loadmaterial('metalAl'), samplekeV);

tube = loaddata('E:\matlab\CTsimulation\physics\tube\tube_spectrumdata_60-150_v1.0.corr');
s = find(tube.KVtag==120);
tube.main = reshape(tube.main, tube.Nsample, []);
spec0 = interp1(tube.main(:, 1), tube.main(:, s*2), samplekeV(:));
% norm
spec0 = spec0./sum(spec0.*samplekeV);

spec1 = exp(-material_Al.mu_total.*2).*spec0;

% constant
electric_charge = 1.602e-19;
epeffect = 0.01;
W = 120*100;
PEscale = 1e-3*W/electric_charge/1000*epeffect;
Eeff = sqrt(sum((spec0.*(samplekeV.^2)))./sum(spec0));

t = 0:0.5:40;
Ptt = exp(-material_Si.mu_total*t).*spec1(:);

t1 = 1.0;
t2 = 39.0;
P1 = exp(-material_Si.mu_total.*t1).*spec1(:);
P2 = exp(-material_Si.mu_total.*t2).*spec1(:);

Ntag = 3;
T = zeros(Ntag-1,1);
for ii = 1:Ntag-1
    Pii = sum(P1)*(1-ii/Ntag) + sum(P2)*(ii/Ntag);
    T(ii) = fzero(@(x) sum(exp(-material_Si.mu_total.*x).*spec1(:))-Pii, 10);
end

T3 = T;

Ntag = 4;
T = zeros(Ntag-1,1);
for ii = 1:Ntag-1
    Pii = sum(P1)*(1-ii/Ntag) + sum(P2)*(ii/Ntag);
    T(ii) = fzero(@(x) sum(exp(-material_Si.mu_total.*x).*spec1(:))-Pii, 10);
end
T4 = T;

Ntag = 6;
T = zeros(Ntag-1,1);
for ii = 1:Ntag-1
    Pii = sum(P1)*(1-ii/Ntag) + sum(P2)*(ii/Ntag);
    T(ii) = fzero(@(x) sum(exp(-material_Si.mu_total.*x).*spec1(:))-Pii, 10);
end
T6 = T;