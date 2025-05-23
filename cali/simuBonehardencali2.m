function bonehardencorr = simuBonehardencali2(SYS, beamhardencorr, corrversion)
% bone harden calibration, slice independent for kneeling effect,
% bonehardencorr = simuBonehardencali2(SYS, beamhardencorr, corrversion)

% default version
if nargin < 3
    corrversion =  'v1.0';
end

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = mergedetector(SYS.detector);
% I know the beamhardencorr has been merged

% paramters
% samplekeV = SYS.world.samplekeV;
focalpos = SYS.source.focalposition(1, :);
detpos = double(detector.position);
Npixel = detector.Npixel;
Nslice = detector.Nslice;
% multi-KV
Nw = SYS.source.Wnumber;

% to cell if not
if ~iscell(beamhardencorr)
    beamhardencorr = num2cell(beamhardencorr);
end

% samplekeV & detector response
det_response = mean(reshape(detector.response, Npixel, []), 1);
det_response = reshape(det_response, Nslice, []);
if strcmpi(SYS.simulation.spectrum, 'Single') && length(det_response)>1
    samplekeV = referencekeV;
    tmp = zeros(Nslice, 1);
    for islice = 1:Nslice
        tmp(islice) =  interp1(world_samplekeV, det_response(islice, :), referencekeV);
    end
    det_response = tmp;
else
    samplekeV = SYS.world.samplekeV;
end
NkeV = length(samplekeV);

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
Dfmu = reshape(Dfmu, Npixel, Nslice, NkeV);
% mean the slices, no
% Dfmu = squeeze(mean(reshape(Dfmu, Npixel, Nslice, []), 2));


% loop source position
for iw = 1:Nw
    % ini
    % m: water, n: bone, q:eff-filter
    % NOTE: the m shall >= beamhardencorr.order
    m = 4;    n = 3;    q = 3;
    m = max(m, beamhardencorr{iw}.order);
    beamhardencorr{iw}.curvematrix = reshape(beamhardencorr{iw}.curvematrix, Nslice, []);
    curvematrix = zeros(Nslice, m*n*q);
    effbeamfilter = zeros(Npixel, Nslice);
    % loop slice
    for islice = 1:Nslice
        detresponse = detspect{iw}(islice, :);
        Dfmu_islice = squeeze(Dfmu(:, islice, :));
        % refernce mu
        mu_wref = interp1(samplekeV, mu_water, beamhardencorr{iw}.referencekeV);
        mu_bref = interp1(samplekeV, mu_bone, beamhardencorr{iw}.referencekeV);
        % response
        Dempty = -log(samplekeV*detresponse')./mu_wref;
        Dfilter = -log(exp(-Dfmu_islice)*(samplekeV.*detresponse)')./mu_wref;
        Deff = Dfilter-Dempty;
        
        % simplify the samples
        [Deff_res, sort_eff] = sort(Deff);
        Nsmp = 20;
        d_effsamp = floor((Deff_res-min(Deff_res))./(max(Deff_res)-min(Deff_res)).*Nsmp);
        [~, i_samp] = unique(d_effsamp);
        Nres = length(i_samp);
        % resample to Nres samples
        Deff_res = Deff_res(i_samp);
        Dfmu_res = Dfmu_islice(sort_eff(i_samp), :);
        
        % projection of air and water+bone
        Pair = -log(exp(-Dfmu_res)*(samplekeV(:).*detresponse))./mu_wref;
        Pres = -log(exp(-repmat(Dwmu, Ndb*Nres, 1) - repmat(repelem(Dbmu, Ndw, 1), Nres, 1) - ...
            repelem(Dfmu_res, Ndw*Ndb, 1))*(samplekeV(:).*detresponse))./mu_wref - repelem(Pair, Ndw*Ndb, 1);
        % Pres = reshape(Pres, Ndw, Ndb, Nres);
        
        % load beamharden curve (surface)
        a = double(beamhardencorr{iw}.curvescale(1));
        b = double(beamhardencorr{iw}.curvescale(2));
        BHmatrix = double(reshape(beamhardencorr{iw}.curvematrix(islice, :), beamhardencorr{iw}.order, []));
        
        D = polyval2dm(BHmatrix, Pres./a, repelem(Deff_res./b, Ndw*Ndb, 1)).*Pres;
        D = reshape(D, Ndw, Ndb, Nres);
        % I know 1st Dbone is 0
        Dw = repmat(D(:,1,:), 1, Ndb, 1);
        Db = repmat(Dbone, Ndw, 1, Nres);
        Df = repmat(reshape(Deff_res, 1, 1, []), Ndw, Ndb, 1);
        
        if islice == 1
            maxD = max(D(:));
            maxDb = max(Db(:));
            maxDf = max(Df(:));
            % will be save in table
            curvescale = [maxD, maxDb, maxDf];
        end
        
        % m: water, n: bone, q:eff-filter
        % X: water+bone, Y: bone, Z: eff-filter
        X = D(:, 2:end, :)./maxD;
        Y = (D(:, 2:end, :)-Dw(:, 2:end, :))./maxDb;
        Z = Df(:, 2:end, :)./maxDf;
        % bone correction target
        Tar = Db(:, 2:end, :)./(D(:, 2:end, :)-Dw(:, 2:end, :));
        % lsqnonlin
        x0 = zeros(1, m*n*q);
        x0(1) = 1;
        options = optimoptions('lsqnonlin', 'Display', 'off');
        x = lsqnonlin(@(x) polyval3dm(reshape(x, m, n, q), X(:), Y(:), Z(:)) - Tar(:), x0, [], [], options);
        
        % will be save in table
        curvematrix(islice, :) = x;
        effbeamfilter(:, islice) = Deff;
    end
    % beamposition
    detpos = reshape(detpos, Npixel, Nslice, []);
    xx=detpos(:,1,1)-focalpos(1);
    yy=detpos(:,1,2)-focalpos(2);
    XYangle = atan2(yy, xx) - pi/2;
    beampos = detector.SID.*sin(XYangle);

    % to table
    bonehardencorr{iw}.ID = corrprm.ID;
    bonehardencorr{iw}.Npixel = Npixel;
    bonehardencorr{iw}.Nslice = Nslice;
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
    bonehardencorr{iw}.effbeamfilter = effbeamfilter;
    % to be manually set
    bonehardencorr{iw}.bonecurvelength = 0;
    bonehardencorr{iw}.bonecurve = [];
end

end
