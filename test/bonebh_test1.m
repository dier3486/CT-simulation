addpath(genpath(pwd));

% % configure files for PX
% calioutputpath = 'E:\matlab\Data\Calibration\BoneBH\';
% system_cfgfile = 'E:\matlab\CT\SINO\PX\system_configure_PXbasic.xml';
% protocol_cfgfile = 'E:\matlab\CT\SINO\PX\simulation\protocol_simulation.xml';
calioutputpath = 'E:\matlab\Data\Calibration\BoneBH\';
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_BoneBH_cali.xml';
protocol_cfgfile = 'E:\matlab\CT\SINO\TM\protocol_beamharden.xml';

% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);

% system configure
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% only play series 3
i_series = 3;
SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
SYS.protocol.seriesindex = i_series;
SYS = loadprotocol(SYS);

% only play 120KV
iw = 3;

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
% 
Nseries = configure.protocol.seriesnumber;
Npixel = SYS.detector.Npixel;
Nslice = max(SYS.detector.slicemerge);      % slice number after merge
Nps = Npixel*Nslice;
KV = SYS.protocol.KV(iw);
samplekeV = SYS.world.samplekeV;
referencekeV = SYS.world.referencekeV;
focalpos = mean(SYS.source.focalposition, 1);
detpos = double(detector.position);
%
det_response = mean(detector.response, 1);
sourcespect = SYS.source.spectrum;
detspect = sourcespect{iw}(:).*det_response(:);

% BH
% bhorder = 3;
% BHcorr = simuBHcali(SYS, bhorder);
% BHcorr = BHcorr{iw};
% or load BH
BHcorr = loaddata('E:\matlab\Data\Calibration\BoneBH\beamharden_Axial_Head_All_120KV300mA_SmallFocalQFS_1SecpRot_v1.0.corr');
BHcorr.curvematrix = reshape(BHcorr.curvematrix, BHcorr.order, []);

% define the bonemodel
bonemodel.material = 'HA';
bonemodel = materialconfigure(bonemodel, samplekeV);

% water sample
Dwater = [2:2:200 204:4:600];
mu_water = SYS.world.water.material.mu_total;
mu_wref = interp1(samplekeV, mu_water, referencekeV);
Dwmu = Dwater(:)*mu_water(:)';
Ndw = length(Dwater);

% bone sample
Dbone = [0:0.5:100 102:2:200];
mu_bone = bonemodel.material.mu_total;
mu_bref = interp1(samplekeV, mu_bone, referencekeV);
Dbmu = Dbone(:)*mu_bone(:)';
Ndb = length(Dbone);

% bowtie and filter
[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
% mean the slices
Dfmu = squeeze(mean(reshape(Dfmu, Npixel, Nslice, []), 2));

Dempty = -log(sum(samplekeV(:).*detspect))./mu_wref;
Dfilter = -log(exp(-Dfmu)*(samplekeV(:).*detspect))./mu_wref;
Deff = Dfilter-Dempty;

% simplify the samples
[Deff_res, sort_eff] = sort(Deff);
Nsmp = 20;
d_effsamp = floor((Deff_res-min(Deff_res))./(max(Deff_res)-min(Deff_res)).*Nsmp);
[~, i_samp] = unique(d_effsamp);
Nres = length(i_samp);
% resample to Nres samples
% Deff_res = repmat(Deff_res(i_samp), 1, Ndw);
% minDeff = min(Deff_res);
% Deff_res = Deff_res-minDeff;
Deff_res = Deff_res(i_samp);
Dfmu_res = Dfmu(sort_eff(i_samp), :);
% Dfilter_res = Dfilter(sort_eff(i_samp));



Pair = -log(exp(-Dfmu_res)*(samplekeV(:).*detspect))./mu_wref;

Pres = -log(exp(-repmat(Dwmu, Ndb*Nres, 1) - repmat(repelem(Dbmu, Ndw, 1), Nres, 1) - ...
       repelem(Dfmu_res, Ndw*Ndb, 1))*(samplekeV(:).*detspect))./mu_wref - repelem(Pair, Ndw*Ndb, 1);
% Pres = reshape(Pres, Ndw, Ndb, Nres);

a = double(BHcorr.curvescale(1));
b = double(BHcorr.curvescale(2));
BHmatrix = double(BHcorr.curvematrix);

D = polyval2dm(BHmatrix, Pres./a, repelem(Deff_res./b, Ndw*Ndb, 1)).*Pres;
D = reshape(D, Ndw, Ndb, Nres);
% I know 1st Dbone is 0
Dw = repmat(D(:,1,:), 1, Ndb, 1);
Db = repmat(Dbone, Ndw, 1, Nres);
Df = repmat(reshape(Deff_res, 1, 1, []), Ndw, Ndb, 1);

maxD = max(D(:));
maxDb = max(Db(:));
maxDf = max(Df(:));

% m: water, n: bone , q:eff-filter
% m=BHcorr.order; n=3; q=3;
% I know m>=BHcorr.order, shall?
m=4; n=3; q=3;

X = D(:, 2:end, :)./maxD;
Y = (D(:, 2:end, :)-Dw(:, 2:end, :))./maxDb;
Z = Df(:, 2:end, :)./maxDf;

Tar = Db(:, 2:end, :)./(D(:, 2:end, :)-Dw(:, 2:end, :));

x0 = zeros(m*n*q, 1);
x0(1) = 1;
options = optimoptions('lsqnonlin','Display','iter');
x = lsqnonlin(@(x) polyval3dm(reshape(x, m, n, q), X(:), Y(:), Z(:)) - Tar(:), x0, [], [], options);

% will be save in table
curvescale = [maxD, maxDb, maxDf];
curvematrix = x;

%
detpos = reshape(detpos, Npixel, Nslice, []);
xx=detpos(:,1,1)-focalpos(1);
yy=detpos(:,1,2)-focalpos(2);
XYangle = atan2(yy, xx) - pi/2;
beampos = detector.SID.*sin(XYangle);
% corrmain = [beampos(:) Deff(:)];

corrprm = parameterforcorr(SYS, 'v1.0');

Bonehardencorr = struct();
Bonehardencorr.ID = corrprm.ID;
Bonehardencorr.Npixel = Npixel;
Bonehardencorr.Nslice = 1;
Bonehardencorr.startslice = corrprm.startslice;
Bonehardencorr.endslice = corrprm.endslice;
Bonehardencorr.slicemerge = corrprm.slicemerge;
Bonehardencorr.focalspot = corrprm.focalspot;
Bonehardencorr.KV = corrprm.KV{iw};
Bonehardencorr.mA = corrprm.mA{iw};
Bonehardencorr.bowtie = corrprm.bowtie;
Bonehardencorr.referencekeV = referencekeV;
Bonehardencorr.refrencemu = mu_wref;
Bonehardencorr.refrencebonemu = mu_bref;
Bonehardencorr.order = [m n q];
Bonehardencorr.matrixsize = m*n*q;
Bonehardencorr.curvescale = curvescale;
Bonehardencorr.curvematrix = curvematrix;
Bonehardencorr.beamposition = beampos;
Bonehardencorr.effbeamfilter = Deff;

bonehardenfile = fullfile(calioutputpath, 'boneharden_test.mat');
save(bonehardenfile, '-struct', 'Bonehardencorr');