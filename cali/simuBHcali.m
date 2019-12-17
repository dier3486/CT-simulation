function BHcorr = simuBHcali(SYS, polyorder, corrversion)
% simulation of bean harden calibration
% BHcorr = simuBHcorr(SYS, response, polyorder)

% default ploly order
if nargin < 2
    polyorder = 4;
end
% default version
if nargin < 4
    corrversion =  'v1.0';
end

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
% paramters
samplekeV = SYS.world.samplekeV;
focalpos = mean(SYS.source.focalposition, 1);
Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;
detpos = double(SYS.detector.position);
% Nsample = length(samplekeV(:));
refrencekeV = SYS.world.refrencekeV;
Nw = SYS.source.Wnumber;

% % test response (debug)
% xr = [486.2832  848.6392  927.3556  734.1626  711.6979];
% spectrange = [20, 150];
% t = linspace(spectrange(1), spectrange(2), length(xr)+1);
% xt = [0; 0; xr(:); 0];
% cs = spline(t, xt);
% response = ppval(cs, samplekeV(:));
% response(samplekeV(:)<20) = 0;
% polyorder = 4;

% spectrums normalize
sourcespect = SYS.source.spectrum;
for iw = 1:Nw
    sourcespect{iw} = double(sourcespect{iw})./sum(double(sourcespect{iw}).*samplekeV);
end
% detector response
detspect = cell(1, Nw);
for iw = 1:Nw
    detspect{iw} = sourcespect{iw}.*detector.spectresponse;
end

% water sample
Dwater = 2:2:600;
mu_water = SYS.world.water.material.mu_total;
mu_wref = interp1(samplekeV, mu_water, refrencekeV);
Dwmu = -Dwater(:)*mu_water(:)';
Ndw = length(Dwater);
% bowtie and filter
[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);

% initial BHcorr
BHcorr = cell(1, Nw);
% parameters for corr
corrprm = parameterforcorr(SYS, corrversion);

% loop Nw
for iw = 1:Nw
    detresponse = detspect{iw}(:);
    
    Dempty = -log(sum(samplekeV(:).*detresponse))./mu_wref;
    Dfilter = -log(exp(-Dfmu)*(samplekeV(:).*detresponse))./mu_wref;
    Deff = Dfilter-Dempty;
    
    % simplify the samples
    [Deff_res, sort_eff] = sort(Deff);
    Nsmp = 200;
    d_effsamp = floor((Deff_res-min(Deff_res))./(max(Deff_res)-min(Deff_res)).*Nsmp);
    [~, i_samp] = unique(d_effsamp);
    Nres = length(i_samp);
    % resample to Nres samples
    Deff_res = repmat(Deff_res(i_samp), 1, Ndw);
    Dfmu_res = Dfmu(sort_eff(i_samp), :);
    Dfilter_res = Dfilter(sort_eff(i_samp));
    
    Pres = exp(-repmat(Dfmu_res, Ndw, 1)+repelem(Dwmu, Nres, 1))*(samplekeV(:).*detresponse);
    Dres = reshape(-log(Pres)./mu_wref, Nres, Ndw) - Dfilter_res;
    Dtarget = Dwater./Dres;
    % fit 2D polynomial
    m = polyorder; n = 4;
    % a = 400; b = 100;
    a = max(Dres(:));
    b = max(Deff_res(:));
    
    x0 = zeros(m*n, 1);
    x0(1) = 1;
    options = optimoptions('lsqnonlin','Display','off');
    x = lsqnonlin(@(x) polyval2dm(reshape(x, m, n), Dres(:)./a, Deff_res(:)./b) - Dtarget(:), x0, [], [], options);
    x = reshape(x, m, n);
    
    % check error (debug)
    % err_res = reshape(polyval2dm(x, Dres(:)./a, Deff_res(:)./b) - Dtarget(:), Nres, Ndw);
    
    % trans x to the polymials of each pixel
    Deff_ply = ones(Npixel*Nslice, n);
    for ii = 2:n
        Deff_ply(:, ii) = (Deff./b).^(ii-1);
    end
    bhpoly = zeros(Npixel*Nslice, m);
    bhpoly(:, m) = Deff_ply*x(1, :)';
    for ii = 1:m-1
        bhpoly(:, m-ii) = Deff_ply*x(ii+1, :)'./a^ii;
    end
    bhpoly(:, 1:m-1) = bhpoly(:, 1:m-1)./bhpoly(:, 2:m);
    
    % normed by mu_weff/log(2)
    bhpoly = bhpoly.*(log(2)/mu_wref);
    
    % double check (debug)
    Dchk = 200;
    Dchk_bh = -log2(exp(-Dfmu - Dchk.*mu_water(:)')*(samplekeV(:).*detresponse)) - Dfilter.*(mu_wref/log(2));
    Dchk_corr = ones(size(Dchk_bh));
    for ii = 1:m
        Dchk_corr = Dchk_corr.*Dchk_bh.*bhpoly(:,ii) + 1.0;
    end
    Dchk_corr = Dchk_corr - 1.0;
    corr_err = (Dchk_corr-Dchk)./Dchk;
    

    
    % slice merge
    [bhpoly, Nmergedslice] = detectorslicemerge(bhpoly, detector, 'mean');
    
    % to table
    BHcorr{iw}.ID = corrprm.ID;
    BHcorr{iw}.Npixel = Npixel;
    BHcorr{iw}.Nslice = Nmergedslice;
    BHcorr{iw}.startslice = corrprm.startslice;
    BHcorr{iw}.endslice = corrprm.endslice;
    BHcorr{iw}.slicemerge = corrprm.slicemerge;
    BHcorr{iw}.focalspot = corrprm.focalspot;
    BHcorr{iw}.KV = corrprm.KV{iw};
    BHcorr{iw}.mA = corrprm.mA{iw};
    BHcorr{iw}.bowtie = corrprm.bowtie;
    BHcorr{iw}.refrencekeV = refrencekeV;
    BHcorr{iw}.refrencemu = mu_wref;
    BHcorr{iw}.order = polyorder;
    BHcorr{iw}.mainsize = Npixel*Nmergedslice*polyorder;
    BHcorr{iw}.main = bhpoly;
end

end