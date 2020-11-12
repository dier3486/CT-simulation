function BHcorr = simuBHcali2(SYS, polyorder, corrversion)
% simulation of bean harden calibration, with kneeling effect, to be merged
% BHcorr = simuBHcorr2(SYS, polyorder)

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
world_samplekeV = SYS.world.samplekeV(:)';
focalpos = mean(SYS.source.focalposition, 1);
Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;
detpos = double(SYS.detector.position);
% Nsample = length(samplekeV(:));
referencekeV = SYS.world.referencekeV;
Nw = SYS.source.Wnumber;

% samplekeV & detector response
det_response = mean(reshape(detector.response, Npixel, []), 1);
det_response = reshape(det_response, Nslice, []);

if strcmpi(SYS.simulation.spectrum, 'Single') && length(det_response)>1
    samplekeV = referencekeV;
    tmp = zeros(Nslice, 1);
    for islice = 1:Nslice
        tmp(islice) =  interp1(world_samplekeV, detector.response(islice, :), referencekeV);
    end
    det_response = tmp;
else
    samplekeV = SYS.world.samplekeV;
end

% spectrums normalize
sourcespect = SYS.source.spectrum;
for iw = 1:Nw
    sourcespect{iw} = double(sourcespect{iw})./sum(double(sourcespect{iw}).*samplekeV);
end
% detector response
detspect = cell(1, Nw);
for iw = 1:Nw
    detspect{iw} = sourcespect{iw}.*det_response;
end

% water sample
Dwater = 2:2:600;
mu_water = SYS.world.water.material.mu_total;
mu_wref = interp1(world_samplekeV, mu_water, referencekeV);
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_water = mu_wref;
end

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
    detresponse = detspect{iw};
    
    Dempty = -log(samplekeV*detresponse')./mu_wref;
    Dfilter = -log(exp(-Dfmu)*(samplekeV.*detresponse)')./mu_wref;
    Deff = Dfilter-Dempty;
    % #1 TBC
    
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
    
    % will be save in table
    curvescale = [a, b];
    curvematrix = x;
    
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
    bhpoly(isnan(bhpoly)) = 0;
    
    % normed by mu_weff/log(2)
    bhpoly = bhpoly.*(log(2)/mu_wref);
    
%     % double check (debug)
%     Dchk = 200;
%     Dchk_bh = -log2(exp(-Dfmu - Dchk.*mu_water(:)')*(samplekeV(:).*detresponse)) - Dfilter.*(mu_wref/log(2));
%     Dchk_corr = ones(size(Dchk_bh));
%     for ii = 1:m
%         Dchk_corr = Dchk_corr.*Dchk_bh.*bhpoly(:,ii) + 1.0;
%     end
%     Dchk_corr = Dchk_corr - 1.0;
%     corr_err = (Dchk_corr-Dchk)./Dchk;
    
    % slice merge
    [bhpoly, Nmergedslice] = detectorslicemerge(bhpoly, detector.Npixel, detector.Nslice, detector.slicemerge, 'mean');
    
    % air rate
    airrate = Deff.*mu_wref./log(2);
    [airrate, ~] = detectorslicemerge(airrate, detector.Npixel, detector.Nslice, detector.slicemerge, 'mean');
    
    % to table
    BHcorr{iw}.ID = corrprm.ID;
    BHcorr{iw}.Npixel = Npixel;
    BHcorr{iw}.Nslice = Nmergedslice;
    BHcorr{iw}.startslice = corrprm.startslice;
    BHcorr{iw}.endslice = corrprm.endslice;
    BHcorr{iw}.mergescale = corrprm.mergescale;
    BHcorr{iw}.focalspot = corrprm.focalspot;
    BHcorr{iw}.KV = corrprm.KV{iw};
    BHcorr{iw}.mA = corrprm.mA{iw};
    BHcorr{iw}.bowtie = corrprm.bowtie;
    BHcorr{iw}.referencekeV = referencekeV;
    BHcorr{iw}.refrencemu = mu_wref;
    BHcorr{iw}.order = polyorder;
    BHcorr{iw}.mainsize = Npixel*Nmergedslice*polyorder;
    BHcorr{iw}.ratesize = Npixel*Nmergedslice;
    BHcorr{iw}.curvescale = curvescale;
    BHcorr{iw}.curvematrix = curvematrix;
    BHcorr{iw}.main = bhpoly;
    BHcorr{iw}.airrate = airrate;
end

end