function bonehardencorr = simuBonehardencali(SYS, beamhardencorr, corrversion)
% bone harden calibration
% bonehardencorr = simuBonehardencali(SYS, beamhardencorr, corrversion)

% default version
if nargin < 3
    corrversion =  'v1.0';
end

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
% paramters
% samplekeV = SYS.world.samplekeV;
focalpos = mean(SYS.source.focalposition, 1);
detpos = double(SYS.detector.position);
Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;
% multi-KV
Nw = SYS.source.Wnumber;

% to cell if not
if ~iscell(beamhardencorr)
    beamhardencorr = {beamhardencorr};
end

% samplekeV & detector response
det_response = mean(detector.response, 1);
if strcmpi(SYS.simulation.spectrum, 'Single') && length(det_response)>1
    samplekeV = SYS.world.referencekeV;
    det_response =  interp1(samplekeV, detector.response, SYS.world.referencekeV);
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

% define the bonemodel
if isfield(SYS.world, 'bonemodel')
    bonemodel = SYS.world.bonemodel;
else
    % default
    bonemodel.material = 'CorticalBone';
    bonemodel = materialconfigure(bonemodel, samplekeV);
end

% water sample (mm)
Dwater = [2:2:200 204:4:600];
mu_water = SYS.world.water.material.mu_total;
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_water = interp1(samplekeV, mu_water, SYS.world.referencekeV);
end
Dwmu = Dwater(:)*mu_water(:)';
Ndw = length(Dwater);

% bone sample (mm)
Dbone = [0:0.5:100 102:2:200];
mu_bone = bonemodel.material.mu_total;
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_bone = interp1(samplekeV, mu_bone, SYS.world.referencekeV);
end
Dbmu = Dbone(:)*mu_bone(:)';
Ndb = length(Dbone);

corrprm = parameterforcorr(SYS, corrversion);
% initial BHcorr
bonehardencorr = cell(1, Nw);

% bowtie and filter
[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
% mean the slices
Dfmu = squeeze(mean(reshape(Dfmu, Npixel, Nslice, []), 2));

% loop source position
for iw = 1:Nw
    % refernce mu
    mu_wref = interp1(samplekeV, mu_water, beamhardencorr{iw}.referencekeV);
    mu_bref = interp1(samplekeV, mu_bone, beamhardencorr{iw}.referencekeV);
    % response
    detresponse = detspect{iw}(:);     
    Dempty = -log(sum(samplekeV(:).*detresponse))./mu_wref;
    Dfilter = -log(exp(-Dfmu)*(samplekeV(:).*detresponse))./mu_wref;
    Deff = Dfilter-Dempty;
    
    % simplify the samples
    [Deff_res, sort_eff] = sort(Deff);
    Nsmp = 20;
    d_effsamp = floor((Deff_res-min(Deff_res))./(max(Deff_res)-min(Deff_res)).*Nsmp);
    [~, i_samp] = unique(d_effsamp);
    Nres = length(i_samp);
    % resample to Nres samples
    Deff_res = Deff_res(i_samp);
    Dfmu_res = Dfmu(sort_eff(i_samp), :);
    
    % projection of air and water+bone
    Pair = -log(exp(-Dfmu_res)*(samplekeV(:).*detresponse))./mu_wref;
    Pres = -log(exp(-repmat(Dwmu, Ndb*Nres, 1) - repmat(repelem(Dbmu, Ndw, 1), Nres, 1) - ...
           repelem(Dfmu_res, Ndw*Ndb, 1))*(samplekeV(:).*detresponse))./mu_wref - repelem(Pair, Ndw*Ndb, 1);
    % Pres = reshape(Pres, Ndw, Ndb, Nres);
    
    % load beamharden curve (surface)
    a = double(beamhardencorr{iw}.curvescale(1));
    b = double(beamhardencorr{iw}.curvescale(2));
    BHmatrix = double(reshape(beamhardencorr{iw}.curvematrix, beamhardencorr{iw}.order, []));

    D = polyval2dm(BHmatrix, Pres./a, repelem(Deff_res./b, Ndw*Ndb, 1)).*Pres;
    D = reshape(D, Ndw, Ndb, Nres);
    % I know 1st Dbone is 0
    Dw = repmat(D(:,1,:), 1, Ndb, 1);
    Db = repmat(Dbone, Ndw, 1, Nres);
    Df = repmat(reshape(Deff_res, 1, 1, []), Ndw, Ndb, 1);

    maxD = max(D(:));
    maxDb = max(Db(:));
    maxDf = max(Df(:));

    % m: water, n: bone, q:eff-filter
    % NOTE: the m shall >= beamhardencorr.order
    m = 4;    n = 3;    q = 3;
    m = max(m, beamhardencorr{iw}.order);
    % X: water+bone, Y: bone, Z: eff-filter
    X = D(:, 2:end, :)./maxD;
    Y = (D(:, 2:end, :)-Dw(:, 2:end, :))./maxDb;
    Z = Df(:, 2:end, :)./maxDf;
    % bone correction target
    Tar = Db(:, 2:end, :)./(D(:, 2:end, :)-Dw(:, 2:end, :));
    % lsqnonlin
    x0 = zeros(m*n*q, 1);
    x0(1) = 1;
    options = optimoptions('lsqnonlin', 'Display', 'off');
    x = lsqnonlin(@(x) polyval3dm(reshape(x, m, n, q), X(:), Y(:), Z(:)) - Tar(:), x0, [], [], options);
    
    % will be save in table
    curvescale = [maxD, maxDb, maxDf];
    curvematrix = x;

    % beamposition
    detpos = reshape(detpos, Npixel, Nslice, []);
    xx=detpos(:,1,1)-focalpos(1);
    yy=detpos(:,1,2)-focalpos(2);
    XYangle = atan2(yy, xx) - pi/2;
    beampos = detector.SID.*sin(XYangle);

    % to table
    bonehardencorr{iw}.ID = corrprm.ID;
    bonehardencorr{iw}.Npixel = Npixel;
    bonehardencorr{iw}.Nslice = 1;
    bonehardencorr{iw}.startslice = corrprm.startslice;
    bonehardencorr{iw}.endslice = corrprm.endslice;
    bonehardencorr{iw}.mergescale = corrprm.mergescale;
    bonehardencorr{iw}.focalspot = corrprm.focalspot;
    bonehardencorr{iw}.KV = corrprm.KV{iw};
    bonehardencorr{iw}.mA = corrprm.mA{iw};
    bonehardencorr{iw}.bowtie = corrprm.bowtie;
    bonehardencorr{iw}.referencekeV = beamhardencorr{iw}.referencekeV;
    bonehardencorr{iw}.refrencemu = mu_wref;
    bonehardencorr{iw}.refrencebonemu = mu_bref;
    bonehardencorr{iw}.order = [m n q];
    bonehardencorr{iw}.matrixsize = m*n*q;
	bonehardencorr{iw}.curvescale = curvescale;
    bonehardencorr{iw}.curvematrix = curvematrix;
    bonehardencorr{iw}.beamposition = beampos;
    bonehardencorr{iw}.effbeamfilter = Deff;
    % to be manually set
    bonehardencorr{iw}.bonecurvelength = 0;
    bonehardencorr{iw}.bonecurve = [];
end

end
