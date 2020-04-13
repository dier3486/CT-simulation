addpath(genpath(pwd));
% inputs
% % configure files for TM
calioutputpath = 'E:\matlab\Data\Calibration\BoneBH\';
system_cfgfile = 'E:\matlab\CT\SINO\TM\system_configure_TM_BoneBH_cali.xml';
protocol_cfgfile = 'E:\matlab\CT\SINO\TM\protocol_beamharden.xml';

% system configure
configure.system = readcfgfile(system_cfgfile);
configure.protocol = readcfgfile(protocol_cfgfile);
configure = configureclean(configure);
Nseries = configure.protocol.seriesnumber;

% system configure
SYS = systemconfigure(configure.system);
SYS = systemprepare(SYS);

% set SYS.output
SYS.output.corrtable = 'bonebeamharden';
% ini the return (corr file name)
BoneBHcalitable = struct();

%need add code to load dfit from beam hardening correciton table
% BHCoef=load('BHCoef.mat');

% loop the body and head
for i_series = 3
    SYS.protocol = protocolconfigure(configure.protocol.series{i_series});
    bowtie = lower(SYS.protocol.bowtie);
    if strcmpi(bowtie, 'empty')
        % skip the BH correction for empty bowtie
        continue;
    end
    SYS.protocol.seriesindex = i_series;
    % load protocol (to SYS)
    SYS = loadprotocol(SYS);
    
    % paramters to use
    Npixel = SYS.detector.Npixel;
    Nslice = max(SYS.detector.slicemerge);      % slice number after merge
    Nps = Npixel*Nslice;
    KV = SYS.protocol.KV;
    Nw = SYS.source.Wnumber;
    
    % ini
    BoneBHcalitable.(bowtie) = cell(1, Nw);
    BoneBHCorr_mat.(bowtie) = cell(1, Nw);
    
    
    
    % loop KV
    for iw = 3
%         dfit = BHCoef.(bowtie).simufilt{iw};
%         x = BHCoef.(bowtie).x{iw};
%         a = BHCoef.(bowtie).a{iw};
%         b = BHCoef.(bowtie).b{iw};
        %load BH cortbl
        BHcor=loaddata('E:\matlab\Data\Calibration\BoneBH\beamharden_Axial_Head_All_120KV300mA_SmallFocalQFS_1SecpRot_v1.0.corr');
        x=reshape(BHcor.curvematrix,[4,4]);
        a=BHcor.curvescale(1);
        b=BHcor.curvescale(2);

        
%         if isempty(dfit)||isempty(x)||isempty(a)||isempty(b)
%             continue;
%         end
        SYS.protocol.KV = KV(iw);
        % reload protocol
        SYS = loadprotocol(SYS);
        % add effect filter
%         Nfilt = length(SYS.collimation.filter);
%         SYS.collimation.filter{Nfilt+1} = struct();
%         SYS.collimation.filter{Nfilt+1}.effect = true;
%         SYS.collimation.filter{Nfilt+1}.thickness = dfit(:);
%         SYS.collimation.filter{Nfilt+1}.material = SYS.collimation.bowtie{1}.material;
        % merge detector slices
        SYS.detector = mergedetector(SYS.detector);
        % BH cali
        BoneBHcorr = BoneBHcali_surffit(SYS,x,a,b);
        % save table
        corrfile = fullfile(calioutputpath, [SYS.output.files.bonebeamharden{1} '.corr']);
        cfgfile = cfgmatchrule(corrfile, SYS.path.IOstandard);
        corrcfg = readcfgfile(cfgfile);
        packstruct(BoneBHcorr{1}, corrcfg, corrfile);

        %save the table file path
        BoneBHcalitable.(bowtie){iw} = corrfile;
        BoneBHCorr_mat.(bowtie){iw}=BoneBHcorr;
    end
end
% save the Corr path in BoneBHcalitable.mat
output_tmp.BoneBHcalitable = BoneBHcalitable;
output_file = fullfile(calioutputpath, 'BHcalitable.mat');
save(output_file, '-struct', 'output_tmp');

% save Corr coefs in BHtemp.mat
output_file = fullfile(calioutputpath, 'BHtemp.mat');
save(output_file, '-struct', 'BoneBHCorr_mat');

