% function corrpoly = simuBHcorr(SYS, response, n)

% test
xr = [486.2832  848.6392  927.3556  734.1626  711.6979];
spectrange = [20, 150];
t = linspace(spectrange(1), spectrange(2), Nx+1);
xt = [0; 0; xr(:); 0];
cs = spline(t, xt);
response = ppval(cs, samplekeV(:));
response(samplekeV(:)<20) = 0;

bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
samplekeV = SYS.world.samplekeV;
focalpos = SYS.source.focalposition(1,:);
Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;
detpos = SYS.detector.position;
Nsample = length(samplekeV(:));
refrencekeV = SYS.world.refrencekeV;

Nbow = length(bowtie(:));

Dwater = 0:1:500;
mu_water = SYS.world.water.material.mu_total;
mu_weff = interp1(samplekeV, mu_water, refrencekeV);
Pwmu = exp(-Dwater(:)*mu_water(:)');
Ndw = length(Dwater);

[Pfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
Deff = -log(exp(-Pfmu)*(samplekeV(:).*response))+log(samplekeV*response);
Deff = reshape(Deff, Npixel, Nslice);
Pfmu = reshape(Pfmu, Npixel, Nslice, Nsample);

ilsice = 1;
% for islice = 1:Nslice
    [Deff_ii, sort_eff] = sort(Deff(:, islice));
    Pfmu_ii = squeeze(Pfmu(sort_eff, islice, :));
    
    Pii = (repmat(Pfmu_ii, Ndw, 1).*repelem(Pwmu, Npixel, 1))*(samplekeV(:).*response);
% end