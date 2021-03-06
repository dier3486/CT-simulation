configure.system = 'E:\matlab\CT\SINO\PG\system_configure_PGcali.xml';

configure = configureclean(configure);
SYS = systemconfigure(configure.system);
if isfield(configure, 'phantom')
    SYS.phantom = phantomconfigure(configure.phantom);
end
SYS = systemprepare(SYS);
protocol.KV = 120;
protocol.bowtie = 2;
protocol.viewnumber = 1;
protocol.collimator = '32x0.625';
SYS.protocol = protocolconfigure(protocol);
SYS = loadprotocol(SYS);

focalpos = SYS.detector.detector_corr.focalposition;
SID = SYS.detector.detector_corr.SID;
FOV = 500;
theta = linspace(-asin(FOV/2/SID), asin(FOV/2/SID), 100);
Nth = size(theta(:), 1);

samplekeV = SYS.world.samplekeV;
water = SYS.world.water;
spectrum = SYS.source.spectrum{1};
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
Nfilter = length(filter);
Nbowtie = length(bowtie);
Nsample = length(samplekeV);

D_filter = zeros(Nth, Nfilter);
mu_filter = zeros(Nfilter, Nsample);
% Dmu = zeros(Nth, Nsample);
for ifilt = 1:Nfilter
    filter_ii = filter{ifilt};
    D_filter(:, ifilt) = filter_ii.thickness./cos(theta);
    mu_filter(ifilt, :) = interp1(filter_ii.material.samplekeV, filter_ii.material.mu_total, samplekeV);
end

D_bowtie = zeros(Nth, Nbowtie);
mu_bowtie = zeros(Nbowtie, Nsample);
% Dmu = zeros(Nth, Nsample);
for ibow = 1:Nbowtie
    bowtie_ii = bowtie{ibow};
    D_bowtie(:, ibow) = interp1(bowtie_ii.anglesample, double(bowtie_ii.bowtiecurve), theta, 'linear', 'extrap');
    mu_bowtie(ibow, :) = interp1(bowtie_ii.material.samplekeV, bowtie_ii.material.mu_total, samplekeV);
end

Dmu = D_filter*mu_filter + D_bowtie*mu_bowtie;

I1 = exp(-Dmu)*samplekeV(:);

Rw = 120;
d = Rw^2 - (SID.*sin(theta)).^2;
d(d<0) = 0;
mu_water = interp1(water.material.samplekeV, water.material.mu_total, samplekeV);
Dw = (sqrt(d(:)).*2)*mu_water;
Iw = exp(-Dmu-Dw)*samplekeV(:);

Iw0 = min(Iw);
% I know
Almaterial = bowtie{2}.material;
mu_Al = interp1(Almaterial.samplekeV, Almaterial.mu_total, samplekeV);
Dal1 = nan(Nth, 1);
for ii = 1:Nth
    if d(ii)==0
        continue;
    end
    Dal1(ii) = fzero(@(x) exp(-D_filter(ii,:)*mu_filter - Dw(ii,:) - x.*mu_Al)*samplekeV(:) - Iw0*1.05, 90);
end
Dal2 = nan(Nth, 1);
for ii = 1:Nth
    if d(ii)==0
        continue;
    end
    Dal2(ii) = fzero(@(x) exp(-D_filter(ii,:)*mu_filter - x.*mu_Al)*samplekeV(:) - I1(ii), 40);
end

Dedg = 32./cos(theta(:));


Rb = 100;
xx1 = Rb.*tan(theta(:)) - Dal1.*sin(theta(:));
yy1 = -Rb + Dal1.*cos(theta(:));
xx2 = Rb.*tan(theta(:)) - Dedg.*sin(theta(:));
yy2 = -Rb + Dedg.*cos(theta(:));


