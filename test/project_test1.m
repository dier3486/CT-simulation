% a projection test script
clear;

% where am I
mainfile = which('CTsimulation');
if isempty(mainfile)
    addpath(genpath('../'));
    mainfile = which('CTsimulation');
    rootpath = fileparts(mainfile);
else
    rootpath = fileparts(mainfile);
end
rootpath = [rootpath '/'];
addpath(genpath(rootpath));

% % configure sample
% configure.system = systemcfgsample(rootpath);
% configure.phantom = phantomcfgsample();
% configure.protocol = protocolcfgsample();
% 
% % output sample xmls
% root = [];
% root.configure = configure;
% struct2xml(root, [rootpath '\system\mod\sample_configure.xml']);
% root = [];
% root.system = configure.system;
% struct2xml(root, [rootpath '\system\mod\sample_system.xml']);
% root = [];
% root.protocol = configure.protocol;
% struct2xml(root, [rootpath '\system\mod\sample_protocol.xml']);

% load configure file
configure_file = 'D:\matlab\ct\BCT16\configure.xml';
configure = readcfgfile(configure_file);

% clean configure
configure = configureclean(configure);

% % output sample xmls 2
% root = [];
% root.configure = configure;
% struct2xml(root, [rootpath 'system\mod\sample_output_configure.xml']);

% get SYS from system configure
SYS = systemconfigure(configure.system);
% phantom
SYS.phantom = phantomconfigure(configure.phantom);

% simulation prepare (load material)
SYS = systemprepare(SYS);

% loop the series
Nseries = configure.protocol.seriesnumber;
for i_series = 1:Nseries
    % to play i-th series
    % load protocol
    SYS.protocol = configure.protocol.series{i_series};
    SYS.protocol.series_index = i_series;
    SYS = loadprotocol(SYS);
    
    % projection (Axial, single, only for test)
    % focal, views, 
    focalposition = SYS.source.focalposition;
    Nfocal = SYS.source.focalnumber;
    Nview = SYS.protocol.viewperrot;
    % Nview = 4;
    viewangle = linspace(0, pi*2, Nview+1);
    viewangle = viewangle(1:end-1) + SYS.protocol.startangle;
    viewangle = reshape(viewangle, Nfocal, []);
    
    Npixel = SYS.detector.Npixel;
    Nslice = SYS.detector.Nslice;
    Np = Npixel * Nslice;
    
    % sigle energy
    samplekeV = SYS.world.refrencekeV;
    Nsample = 1;
    energynorm = 1/SYS.world.refrencekeV;
    % ini Dmu
    Dmu = zeros(Np*Nview, Nsample);
    % projection on bowtie and filter
    % subfunction: bowtie-projection
    xx = SYS.detector.position(:,1) - focalposition(:,1)';
    yy = SYS.detector.position(:,2) - focalposition(:,2)';
    zz = SYS.detector.position(:,3) - focalposition(:,3)';
    detangle = atan2(yy, xx) - pi/2;
    detZscale = sqrt(yy.^2+zz.^2)./yy;  
    % bowtie(s)
    Nbowtie = length(SYS.collimation.bowtie(:));
    for ibow = 1:Nbowtie
        bowtie = SYS.collimation.bowtie{ibow};
        if isempty(bowtie.bowtiecurve)
            % empty bowtie
            continue;
        end
        % D
        D_bowtie = interp1(bowtie.anglesample, bowtie.bowtiecurve, detangle);
        D_bowtie = D_bowtie.*detZscale;
        % mu
        if Nsample == 1
            % sigle energy
            mu_bowtie = interp1(bowtie.material.samplekeV, bowtie.material.mu_total, samplekeV);
        else
            mu_bowtie = bowtie.material.mu_total;
        end
        % + to Dmu
        Dmu = Dmu + repmat(D_bowtie(:)*mu_bowtie, Nview/Nfocal, 1);
    end
    % filter(s)
    Nfilter = length(SYS.collimation.filter(:));
    Dfscale = (sqrt(xx.^2+yy.^2+zz.^2)./yy);
    for ifil = 1:Nfilter
        filter = SYS.collimation.filter{ifil};
        % D
        D_filter = Dfscale.*filter.thickness;
        % mu 
        if Nsample == 1
            % sigle energy
            mu_filter = interp1(filter.material.samplekeV, filter.material.mu_total, samplekeV);
        else
            mu_filter = filter.material.mu_total;
        end
        % + to Dmu
        Dmu = Dmu + repmat(D_filter(:)*mu_filter, Nview/Nfocal, 1);
    end
    % P(osibility) of air
    P_air = exp(-Dmu(1:Np*Nfocal))*samplekeV';
    % no detector response employed
        
    % projection in phantoms
    % subfunction: phantom-projection
    for iobj = 1:SYS.phantom.Nobject
        parentobj = SYS.phantom.object_tree(iobj);
        object_i = SYS.phantom.object{iobj};
        if Nsample == 1
            % sigle energy
            mu_i = interp1(object_i.material.samplekeV, object_i.material.mu_total, samplekeV);
        else
            mu_i = object_i.material.mu_total;
        end
        if parentobj>0
            if Nsample == 1
                % sigle energy
                mu_parent = interp1(SYS.phantom.object{parentobj}.material.samplekeV, ...
                    SYS.phantom.object{parentobj}.material.mu_total, samplekeV);
            else
                mu_parent = SYS.phantom.object{parentobj}.material.mu_total;
            end
            mu_i = mu_i - mu_parent;
        end
        [D_i, L] = intersection(focalposition, SYS.detector.position, object_i, 'views-ray', viewangle, 0);
        Dmu = Dmu + D_i(:)*mu_i;
    end
    % P(osibility)
    P = exp(-Dmu)*samplekeV';
    P = reshape(P, Np, Nview);
    % no detector response employed
    
    % distance curse
    pixel_area = 1;
    distscale = pixel_area./(L.^2.*(pi*4));
    P_air = P_air.*distscale;
    P = P.*distscale;
    
    % measurement parameters, DCB
    electric_charge = 1.602e-19;
    KV = SYS.protocol.KV;
    mA = SYS.protocol.mA;
    mA_air = mA;
    W = KV*mA;
    W_air = W;
    T = SYS.protocol.rotationspeed/SYS.protocol.viewperrot.*1e6;
    gain = 0.1;
    Z0 = 16384;
    maxanglecode = 69120;
    
    % Intansity
    Pscale = (T*1e-6*W/electric_charge*energynorm/1000).*gain;
    Intensity = P.*Pscale + Z0;
    
    % to put in rawdata struct
    Intensity = num2cell(Intensity, 1);
    readingnumber = num2cell(1:Nview, 1);
    angleencoder = mod(round(viewangle./(pi*2/maxanglecode)), maxanglecode);
    angleencoder = num2cell(angleencoder, 1);
    rawdataversion = [1 0];
    statusflag = hex2dec('8000');
    
    % rawdata struct
    raw = struct();
    raw(Nview) = struct();
    [raw(:).Package_Version] = deal(rawdataversion);
    [raw(:).Status_Flag] = deal(statusflag);
    [raw(:).Reading_Number] = readingnumber{:};
    [raw(:).Angle_encoder] = angleencoder{:};
    [raw(:).Integration_Time] = deal(T*125);
    [raw(:).KV] = deal(KV);
    [raw(:).mA] = deal(mA);
    [raw(:).Start_Slice] = deal(1);
    [raw(:).End_Slice] = deal(Nslice);
    [raw(:).Raw_Data_Size] = deal(Np*3);
    [raw(:).Slice_Number] = deal(Nslice);
    [raw(:).Raw_Data] = Intensity{:};
    
    % rawdata output
    % file name
    rawdatafile = [SYS.output.path SYS.output.files.rawdata '_' SYS.output.rawdataversion '.raw'];
    % find the format configure file
    rawcfgfile = cfgmatchrule(rawdatafile, SYS.path.IOstandard, SYS.output.rawdataversion);
    rawcfg = readcfgfile(rawcfgfile);
    % pack the data
	packstruct(raw, rawcfg, rawdatafile);
    
    % air calibration
    rawair = log2(P_air.*Pscale.*(mA_air/mA)) - log2(T*125);
    
    % air calibration table
    % corr table base
    aircorr_cfgfile = [rootpath, 'IO/standard/air_corr_v1.0.xml'];
    aircorr_basefile = [rootpath, 'IO/standard/air_sample_v1.0.corr'];
    aircorr_base = loadbindata(aircorr_basefile, aircorr_cfgfile);
    % fill up
    aircorr = aircorr_base;
    aircorr.ID = [0 0 1 0];
    aircorr.main = rawair(:);
    aircorr.mainsize = size(rawair(:),1);
    aircorr.Nsection = 1;
    aircorr.reference = 1.0;
    aircorr.Nslice = Nslice;
    aircorr.Npixel = Npixel;
    aircorr.KV = KV;
    aircorr.mA = mA_air;
    
    % output air corr table
	aircorrfile = [SYS.output.path SYS.output.files.aircorr '_v1.0' '.corr'];
    aircfgfile = cfgmatchrule(aircorrfile, SYS.path.IOstandard, 'v1.0');
    aircfg = readcfgfile(aircfgfile);
    packstruct(aircorr, aircfg, aircorrfile);
end


