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
detpos = double(SYS.detector.position);
Nsample = length(samplekeV(:));
refrencekeV = SYS.world.refrencekeV;

Nbow = length(bowtie(:));

Dwater = 2:2:600;
mu_water = SYS.world.water.material.mu_total;
mu_weff = interp1(samplekeV, mu_water, refrencekeV);
Dwmu = -Dwater(:)*mu_water(:)';
Ndw = length(Dwater);

[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
Dempty = -log(samplekeV*response)./mu_weff;
Dfilter = reshape(-log(exp(-Dfmu)*(samplekeV(:).*response))./mu_weff, Npixel, Nslice);
Deff = Dfilter-Dempty;
Dfmu = reshape(Dfmu, Npixel, Nslice, Nsample);

m = 4; n = 4;
a = 400; b = 100;
bhpoly = zeros(Npixel*Nslice, m);

islice = 1;
% for islice = 1:Nslice
    [Deff_ii, sort_eff] = sort(Deff(:, islice));
    Dfmu_ii = squeeze(Dfmu(sort_eff, islice, :));
    % simplify
    Nsmp = 200;
    d_effsamp = floor((Deff_ii-min(Deff_ii))./(max(Deff_ii)-min(Deff_ii)).*Nsmp);
    [d_samp, i_samp] = unique(d_effsamp);
    Nrep = length(i_samp);
    Deff_ii = repmat(Deff_ii(i_samp), 1, Ndw);
    Dfmu_ii = Dfmu_ii(i_samp, :);
    Dfilter_ii = Dfilter(sort_eff, islice);
    Dfilter_ii = Dfilter_ii(i_samp);
      
    Pii = exp(-repmat(Dfmu_ii, Ndw, 1)+repelem(Dwmu, Nrep, 1))*(samplekeV(:).*response);
    Dii = reshape(-log(Pii)./mu_weff, Nrep, Ndw) - Dfilter_ii;
    Dtarget = Dwater./Dii;
    
    x0 = zeros(m*n, 1);
    x0(1) = 1;
    
    x = lsqnonlin(@(x) polyval2dm(reshape(x, m, n), Dii(:)./a, Deff_ii(:)./b) - Dtarget(:), x0);
    x = reshape(x, m, n);
    % check error
    err_ii = reshape(polyval2dm(x, Dii(:)./a, Deff_ii(:)./b) - Dtarget(:), Nrep, Ndw);
    
    Deff_ply = ones(Npixel, n);
    for ii = 2:n
        Deff_ply(:, ii) = (Deff(:, islice)./b).^(ii-1);
    end
    index_ii = (1:Npixel) + (islice-1)*Npixel;
    bhpoly(index_ii, m) = Deff_ply*x(1, :)';
    for ii = 1:m-1
        bhpoly(index_ii, m-ii) = Deff_ply*x(ii+1, :)'./a^ii;
    end
    bhpoly(:, 1:m-1) = bhpoly(:, 1:m-1)./bhpoly(:, 2:m);
    
    % check
    Dchk = 180;
    Pchk = exp(-squeeze(Dfmu(:,islice,:))-Dchk.*mu_water(:)')*(samplekeV(:).*response);
    Dchk_bh = -log(Pchk)./mu_weff - Dfilter(:, islice);
    Dchk_corr = ones(size(Dchk_bh));
    for ii=1:m
        Dchk_corr = Dchk_corr.*Dchk_bh.*bhpoly(index_ii, ii) + 1.0;
    end
    Dchk_corr = Dchk_corr - 1.0;
    
    
% end