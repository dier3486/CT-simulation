function Bonehardencorr = BoneBHcali_3Dfit(SYS,corrversion)

% default version
if nargin < 2
    corrversion =  'v1.0';
end

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
% paramters
samplekeV = SYS.world.samplekeV;
focalpos = mean(SYS.source.focalposition, 1);
detpos = double(SYS.detector.position);

Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;


%后续需要改为使用水硬化校正表存储的mu_ref
referencekeV = SYS.world.referencekeV;
Nw = SYS.source.Wnumber;

% samplekeV & detector response
det_response = mean(detector.response, 1);
if strcmpi(SYS.simulation.spectrum, 'Single') && length(det_response)>1
    samplekeV = referencekeV;
    det_response =  interp1(samplekeV, detector.response, referencekeV);
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

% BH
% bhorder = 3;
% BHcorr = simuBHcali(SYS, bhorder);
% BHcorr = BHcorr{iw};
% or load BH
BHcorr = loaddata('E:\data\simulation\PG\beamharden_Axial_Head_32x0.625_140KV300mA_SmallFocalQFS_1SecpRot_v1.0.corr');
BHcorr.curvematrix = reshape(BHcorr.curvematrix, BHcorr.order, []);

% define the bonemodel
% bonemodel.material = 'HA';
bonemodel.material = 'CorticalBone';
bonemodel = materialconfigure(bonemodel, samplekeV);

% water sample
Dwater = [2:2:200 204:4:600];
mu_water = SYS.world.water.material.mu_total;
mu_wref = interp1(samplekeV, mu_water, referencekeV);
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_water = mu_wref;
end
Dwmu = Dwater(:)*mu_water(:)';
Ndw = length(Dwater);

% bone sample
Dbone = [0:0.5:100 102:2:200];
mu_bone = bonemodel.material.mu_total;
mu_bref = interp1(samplekeV, mu_bone, referencekeV);
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_bone = mu_bref;
end
Dbmu = Dbone(:)*mu_bone(:)';
Ndb = length(Dbone);

corrprm = parameterforcorr(SYS, corrversion);
% initial BHcorr
Bonehardencorr = cell(1, Nw);

% bowtie and filter
[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
% mean the slices
Dfmu = squeeze(mean(reshape(Dfmu, Npixel, Nslice, []), 2));

% [m,n] = size(x);
% 
% %计算中心层，每个通道射线路径距离ISO中心的距离chpos
% xpos=reshape(detpos(:,1),SYS.detector.Npixel,SYS.detector.Nslice);
% xpos=xpos(:,1);
% chPos=xpos*SYS.detector.SID/SYS.detector.SDD;
% 
% %将Dfmu重采样为以mm为单位间隔的数组Dfmu_res
% chPos_min=floor(min(chPos));
% chPos_max=ceil(max(chPos));
% Dfmu=Dfmu(Npixel*Nslice/2+1:Npixel*Nslice/2+Npixel,:);
% Dfmu_res=interp1(chPos,Dfmu,chPos_min:0.5:chPos_max,'spline');
% %如果固定正投影通道数为1000，则Dfmu_res按下面式子计算
% % Dfmu_res=interp1(chPos,Dfmu,-249.5:0.5:250,'spline');
% Nres=length(Dfmu_res);



% loop source position
for iw = 1:Nw
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
    % Deff_res = repmat(Deff_res(i_samp), 1, Ndw);
    % minDeff = min(Deff_res);
    % Deff_res = Deff_res-minDeff;
    Deff_res = Deff_res(i_samp);
    Dfmu_res = Dfmu(sort_eff(i_samp), :);
    % Dfilter_res = Dfilter(sort_eff(i_samp));

    Pair = -log(exp(-Dfmu_res)*(samplekeV(:).*detresponse))./mu_wref;

    Pres = -log(exp(-repmat(Dwmu, Ndb*Nres, 1) - repmat(repelem(Dbmu, Ndw, 1), Nres, 1) - ...
           repelem(Dfmu_res, Ndw*Ndb, 1))*(samplekeV(:).*detresponse))./mu_wref - repelem(Pair, Ndw*Ndb, 1);
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

    % m: bone, n: water, q:eff-filter
    % m=BHcorr.order; n=3; q=3;
    % I know n>=BHcorr.order, shall?
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

%     corrprm = parameterforcorr(SYS, 'v1.0');
% 
%     Bonehardencorr = struct();
    Bonehardencorr{iw}.ID = corrprm.ID;
    Bonehardencorr{iw}.Npixel = Npixel;
    Bonehardencorr{iw}.Nslice = 1;
    Bonehardencorr{iw}.startslice = corrprm.startslice;
    Bonehardencorr{iw}.endslice = corrprm.endslice;
    Bonehardencorr{iw}.mergescale = corrprm.mergescale;
    Bonehardencorr{iw}.focalspot = corrprm.focalspot;
    Bonehardencorr{iw}.KV = corrprm.KV{1};
    Bonehardencorr{iw}.mA = corrprm.mA{1};
    Bonehardencorr{iw}.bowtie = corrprm.bowtie;
    Bonehardencorr{iw}.referencekeV = referencekeV;
    Bonehardencorr{iw}.refrencemu = mu_wref;
    Bonehardencorr{iw}.refrencebonemu = mu_bref;
    Bonehardencorr{iw}.order = [m n q];
    Bonehardencorr{iw}.matrixsize = m*n*q;
    Bonehardencorr{iw}.curvescale = curvescale;
    Bonehardencorr{iw}.curvematrix = curvematrix;
    Bonehardencorr{iw}.beamposition = beampos;
    Bonehardencorr{iw}.effbeamfilter = Deff;
   
    % slice merge
%     [bhpoly, Nmergedslice] = detectorslicemerge(bhpoly, detector.Npixel, detector.Nslice, detector.slicemerge, 'mean');
    % to table

end

end
