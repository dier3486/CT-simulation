% run after expbgback_test1.m

% use the effect filter to get BHcorr

% CT system
configure_file = 'E:\matlab\calibration\system\configure_cali.xml';
% read configure file
configure = readcfgfile(configure_file);
% load configure
configure = configureclean(configure);
% system configure
SYS = systemconfigure(configure.system);
% phantom configure
SYS.phantom = phantomconfigure(configure.phantom);
% simulation prepare (load materials)
SYS = systemprepare(SYS);

% response 
rsp = load('./testdata/response_1219.mat');
% effect filter
efffilt = load('./testdata/bowtiefit_1220.mat');

% output path
outpath = './testdata/';
corrversion = 'v1.0';

% ini return
Cali = struct();
Cali(8) = struct();
% to loop the bowtie and KV
KV = [80 100 120 140];
bowtie = {'Body', 'Head'};
for ibow = 1:length(bowtie)
    for iKV = 1:length(KV)
        SYS.protocol = protocolconfigure(configure.protocol.series{1});
        SYS.protocol.series_index = 1;
        % change bowtie
        SYS.protocol.bowtie = bowtie{ibow};
        % change KV
        SYS.protocol.KV = KV(iKV);
        % load protocol (to SYS)
        SYS = loadprotocol(SYS);
        % add effect filter
        Nfilt = length(SYS.collimation.filter);
        SYS.collimation.filter{Nfilt+1} = struct();
        SYS.collimation.filter{Nfilt+1}.effect = true;
        thickness = efffilt.bowtiefit{iKV}(:,:, ibow);
        SYS.collimation.filter{Nfilt+1}.thickness = thickness(:);
        SYS.collimation.filter{Nfilt+1}.material = SYS.collimation.bowtie{1}.material;
        % add response
        SYS.detector.spectresponse = rsp.response(:)';
        
        % BH cali
        BHcorr = simuBHcali(SYS, 4);
        index_cal = (ibow-1)*length(KV) + iKV;
        Cali(index_cal).BHcorr = BHcorr{1};
        Cali(index_cal).KV = KV(iKV);
        Cali(index_cal).bowtie = bowtie{ibow};
        
        % save table
        tablename = ['beamharden_' bowtie{ibow} '_' num2str(KV(iKV)) '_v1.0.corr'];
        corrfile = fullfile(outpath, tablename);
        cfgfile = cfgmatchrule(corrfile, SYS.path.IOstandard);
        corrcfg = readcfgfile(cfgfile);
        packstruct(Cali(index_cal).BHcorr, corrcfg, corrfile);
    end
end

% the returns are
% Cali
