% script for response fitting

% bladeexperiment

samplekeV = SYS.world.samplekeV;
mu_Cu = SYS.phantom.object{1}.material.mu_total;
mu_ref = interp1(samplekeV, mu_Cu, SYS.world.refrencekeV);

sampixel = 400:460;
samslice = 24:39;
samindex = sampixel(:) + (samslice-1).*Npixel;
samindex = samindex(:);

Dexp = zeros(Nbld, Nw);
for iw = 1:Nw
    Dexp(:, iw) = real(-log(mean(expbld{iw}(samindex, :), 1))' + log(mean(expair{iw}(samindex))))./mu_ref;
end
Dexp = Dexp(:);
Savl = Dexp<10.0;

Abld = zeros(Nbld*Nw, Nsmp);
Aair = zeros(Nbld*Nw, Nsmp);
for iw = 1:Nw
    index  = (1:Nbld) + (iw-1)*Nbld;
    Aair(index, :) = repmat(mean(Data.Pair{iw}(samindex, :), 1), Nbld, 1);
    Abld(index, :) = squeeze(mean(Pbld{iw}(samindex, :, :), 1))';
end

r0 = ones(size(samplekeV));
D0 = (-log(Abld*(r0(:).*samplekeV(:))) + log(Aair*(r0(:).*samplekeV(:))))./mu_ref;


spectrange = [0 , 140];
x0 = [ 1 1 1 1 1 1 1 ];


Abld_fit = Abld(Savl, :).*samplekeV;
Aair_fit = Aair(Savl, :).*samplekeV;
Dexp_fit = Dexp(Savl);

x1 = x0;
x = lsqnonlin(@(x) spectfitfun(x, Abld_fit, Aair_fit, spectrange, samplekeV, mu_ref, Dexp_fit), x1);

xt = [0; abs(x(:));];
Nx = length(xt);
t = linspace(spectrange(1), spectrange(2), Nx);
r2 = pchip(t, xt, samplekeV(:));
% r2(samplekeV(:)<spectrange(1)) = 0;
% r2 = r2.*Cii.GOS;
D2 = (-log(Abld*(r2(:).*samplekeV(:))) + log(Aair*(r2(:).*samplekeV(:))))./mu_ref;

figure;
plot(samplekeV, r2);

figure;
hold on
plot(blades, reshape(D2,[],4));
plot(blades, reshape(Dexp,[],4), 'o');
axis([0 20 0 12]);
axis equal
grid on

