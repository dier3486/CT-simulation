BoneBHCorr=load(strCortblpath);


%找出每个通道需要使用的校正表位置
fChannelSpace=FPStruct.fSysFov / FPStruct.nChaNum;
nChannelNum=FPStruct.nChaNum;
fMidChannel=FPStruct.fCenCha;
fChannelPos = ((1:nChannelNum)-fMidChannel)*fChannelSpace;

% prepare
BoneBHStruct.efffilter = interp1(BoneBHCorr.beamposition, BoneBHCorr.effbeamfilter, fChannelPos, 'linear', 'extrap')./BoneBHCorr.curvescale(3);
BoneBHStruct.efffilter = BoneBHStruct.efffilter(:);
BoneBHStruct.curvematrix = reshape(BoneBHCorr.curvematrix, BoneBHCorr.order);
BoneBHStruct.Dscale = BoneBHCorr.curvescale(1);
BoneBHStruct.bonescale = BoneBHCorr.curvescale(2);
BoneBHStruct.mubonmuw = BoneBHCorr.refrencebonemu/BoneBHCorr.refrencemu;
BoneBHStruct.HCscale = 1000;

% corr

D = ImgProjTemp./BoneBHStruct.HCscale./BoneBHStruct.Dscale;
Db = BoneProjTemp./BoneBHStruct.HCscale./BoneBHStruct.bonescale;

