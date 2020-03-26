function BoneBHcor = BoneBHcali_surffit(SYS,x,a,b,corrversion)

% default version
if nargin < 5
    corrversion =  'v1.0';
end

% system components
bowtie = SYS.collimation.bowtie;
filter = SYS.collimation.filter;
detector = SYS.detector;
% paramters
world_samplekeV = SYS.world.samplekeV;
focalpos = mean(SYS.source.focalposition, 1);
Npixel = SYS.detector.Npixel;
Nslice = SYS.detector.Nslice;
detpos = double(SYS.detector.position);

%后续需要改为使用水硬化校正表存储的mu_ref
referencekeV = SYS.world.referencekeV;
Nw = SYS.source.Wnumber;

% samplekeV & detector response
det_response = mean(detector.response, 1);
if strcmpi(SYS.simulation.spectrum, 'Single') && length(det_response)>1
    samplekeV = referencekeV;
    det_response =  interp1(world_samplekeV, detector.response, referencekeV);
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
Dwater = 0:0.5:600;
mu_water = SYS.world.water.material.mu_total;
mu_wref = interp1(world_samplekeV, mu_water, referencekeV);
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_water = mu_wref;
end

Dwmu = -Dwater(:)*mu_water(:)';
Ndw = length(Dwater);

% bone sample
% Dbone = 0:2:500;
Dbone = 0:0.2:50;
mu_bone = SYS.world.bone.material.mu_total;
mu_bref = interp1(world_samplekeV, mu_bone, referencekeV);
if strcmpi(SYS.simulation.spectrum, 'Single')
    mu_bone = mu_bref;
end

Dbmu = -Dbone(:)*mu_bone(:)';
Ndb = length(Dbone);

corrprm = parameterforcorr(SYS, corrversion);
% initial BHcorr
BoneBHcor = cell(1, Nw);

% bowtie and filter
[Dfmu, ~] = flewoverbowtie(focalpos, detpos, bowtie, filter, samplekeV);
[m,n] = size(x);

%计算中心层，每个通道射线路径距离ISO中心的距离chpos
xpos=reshape(detpos(:,1),864,64);
xpos=xpos(:,1);
chPos=xpos*SYS.detector.SID/SYS.detector.SDD;

%将Dfmu重采样为以mm为单位间隔的数组Dfmu_res
chPos_min=floor(min(chPos));
chPos_max=ceil(max(chPos));
Dfmu=Dfmu(Npixel*Nslice/2+1:Npixel*Nslice/2+Npixel,:);
Dfmu_res=interp1(chPos,Dfmu,chPos_min:0.5:chPos_max,'spline');
%如果固定正投影通道数为1000，则Dfmu_res按下面式子计算
% Dfmu_res=interp1(chPos,Dfmu,-249.5:0.5:250,'spline');
Nres=length(Dfmu_res);

% loop source position
for iw = 1:Nw
    detresponse = detspect{iw}(:); 
    Dempty = -log(sum(world_samplekeV(:).*detresponse))./mu_wref;
    Dfilter = -log(exp(-Dfmu_res)*(samplekeV(:).*detresponse))./mu_wref;
    Deff = Dfilter-Dempty;     
    
    %  beamhardening correciton prepare
    Deff_ply = ones(Nres, n);
    for ii = 2:n
        Deff_ply(:, ii) = (Deff./b).^(ii-1);
    end
    
    %计算不同filter厚度时，水硬化校正曲线系数。
    BHCoef = Deff_ply * x';
  
    % 设置骨硬化校正多项式系数矩阵大小
    u = 3; v = 3;   
    cor=zeros(u,v,Nres);   
    
    %Loop 不同的滤过厚度 Dfmu_res 
    for i=1:Nres
        Pres = exp(-Dfmu_res(i,:)+repelem(Dwmu, Ndb, 1)+repmat(Dbmu, Ndw, 1))*(samplekeV(:).*detresponse);
        Dres = (-log(Pres)./mu_wref - Dfilter(i))/a;
        %使用x对Dres_B(:,i)进行水硬化校正
        Dres_poly=ones(Ndb*Ndw, m);
        for j=1:m
            Dres_poly(:,j) = Dres.^j;
        end
        tempDres = sum(BHCoef(i,:).*Dres_poly,2)*a;
        D = double(reshape(tempDres,[Ndb, Ndw]));

        [Bid,Wid]=ndgrid(Dbone,Dwater);
        
        Dd_w = D-Wid;
        Db = Bid;
        x0 = zeros(u*v, 1);
        x0(1) = 0;
        options = optimoptions('lsqnonlin','Display','off');

        xb = lsqnonlin(@(x) polyval2dm(reshape(x, u, v),Dd_w,D) - Db , x0, [], [], options);
        cor(:,:,i) = reshape(xb, u, v);
        i
    end    
    % slice merge
%     [bhpoly, Nmergedslice] = detectorslicemerge(bhpoly, detector.Npixel, detector.Nslice, detector.slicemerge, 'mean');
    % to table
    BoneBHcor{iw}.ID = corrprm.ID;
    BoneBHcor{iw}.startslice = corrprm.startslice;
    BoneBHcor{iw}.endslice = corrprm.endslice;
    BoneBHcor{iw}.slicemerge = corrprm.slicemerge;
    BoneBHcor{iw}.focalspot = corrprm.focalspot;
    BoneBHcor{iw}.KV = corrprm.KV{iw};
    BoneBHcor{iw}.mA = corrprm.mA{iw};
    BoneBHcor{iw}.bowtie = corrprm.bowtie;    
    BoneBHcor{iw}.Nres = Nres;    
    BoneBHcor{iw}.referencekeV = referencekeV;
    BoneBHcor{iw}.refrencemu = mu_wref;
    BoneBHcor{iw}.order = 2;
    BoneBHcor{iw}.mainsize = Nres*u*v;
    BoneBHcor{iw}.main = cor;

end

end
