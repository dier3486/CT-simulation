% script for response fitting

% bladeexperiment

samplekeV = SYS.world.samplekeV;
mu_Cu = SYS.phantom.object{1}.material.mu_total;
mu_ref = interp1(samplekeV, mu_Cu, SYS.world.referencekeV);

mid_u = 427.75;
% sampixel = 400:460;
% samslice = 24:39;
sampixel = round(mid_u)-30 : round(mid_u)+30;
samslice = 9:24;

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


spectrange = [20 , 140];
x0 = [ 1 1 1 1 1 1 ];


Abld_fit = Abld(Savl, :).*samplekeV;
Aair_fit = Aair(Savl, :).*samplekeV;
Dexp_fit = Dexp(Savl);

% x1 = x0;
x1 = [12.4420   16.8111   16.0373   17.9764   17.1018   13.7785];
x = lsqnonlin(@(x) spectfitfun(x, Abld_fit, Aair_fit, spectrange, samplekeV, mu_ref, Dexp_fit), x1);

xt = [abs(x(:));];
Nx = length(xt);
t = linspace(spectrange(1), spectrange(2), Nx);
r2 = pchip(t, xt, samplekeV)./mean(xt);
% r2(samplekeV(:)<spectrange(1)) = 0;
% r2 = r2.*Cii.GOS;
D2 = (-log(Abld*(r2(:).*samplekeV(:))) + log(Aair*(r2(:).*samplekeV(:))))./mu_ref;

figure;
plot(samplekeV, r2);
hold on
plot(t, xt./mean(xt), 'o');
grid on

figure;
hold on
plot(blades, reshape(D2,[],4));
plot(blades, reshape(Dexp,[],4), 'o');
axis([0 20 0 12]);
axis equal
grid on

% return
resp.samplekeV = samplekeV;
resp.response = r2;

