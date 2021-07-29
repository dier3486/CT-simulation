% test of Compton formula
% run a system configure script first

samplekeV = SYS.world.samplekeV;
f0 = SYS.source.spectrum{1};
% f0(f0==0) = nan;

t = 1/16;

samplekeV2 = samplekeV(1):t:samplekeV(end);

mec2 = 511;

Nth = 1000;
theta = linspace(-pi, pi, Nth)';
ctm = (1-cos(theta))./mec2;

invg = 1./(1./samplekeV2 - ctm);
finvg = interp1(samplekeV, f0, invg);
finvg = fillmissing(finvg, 'const', 0);
dginvg = 1./(ctm.*invg + 1).^2;

f1 = finvg./dginvg;
y1 = f1*samplekeV2(:).*t;

% Oxygen
Z = 8;
Wz = 15.9994;

alphare2 = 7.940775e-30*1e4; % m^2 -> cm^2
Watom = 1.66053904e-27*1e3;  % kg -> g
lambda_A0 = 12398.520;

KeV0 = 50;
sigma0_coh = 2.166E-02;
sigma0_inc = 1.609E-01;
dataO = readmatrix('E:\matlab\CTsimulation\test\O_FS.csv');

q = 1./(1+ctm.*KeV0);
dsigmaKN = (alphare2/2/Watom/Wz).*q.^2.*(q+1./q-sin(theta).^2);
dsigmaR = (alphare2/2/Watom/Wz).*(1+cos(theta).^2);
x_theta = abs(sin(theta./2)).*(KeV0*1e3/lambda_A0);
F_theta = interp1(dataO(:,1), dataO(:,2), x_theta, 'pchip');
S_theta = interp1(dataO(:,1), dataO(:,3), x_theta, 'pchip');

w = [1/2 ones(1, Nth/2-1)].*(pi*2/(Nth-1));
Sinc = dsigmaKN.*S_theta.*abs(sin(theta)).*(pi*2);
Scoh = dsigmaR.*F_theta.^2.*abs(sin(theta)).*(pi*2);

sig_inc = w*Sinc(1:Nth/2);
sig_coh = w*Scoh(1:Nth/2);

figure;
hold on;
plot(dsigmaKN.*S_theta./Z.*cos(theta), dsigmaKN.*S_theta./Z.*sin(theta));
plot(dsigmaR.*F_theta.^2./Z^2.*cos(theta), dsigmaR.*F_theta.^2./Z^2.*sin(theta));
axis equal
