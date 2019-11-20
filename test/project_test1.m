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
addpath(genpath(rootpath));

configure.system = systemcfgsample([rootpath '\']);
configure.phantom = phantomcfgsample();
configure.protocol = protocolcfgsample();

% output sample xmls
% root = [];
% root.configure = configure;
% struct2xml(root, 'D:\matlab\CTsimulation\system\mod\sample_configure.xml');
root = [];
root.system = configure.system;
struct2xml(root, [rootpath 'system\mod\sample_system.xml']);
root = [];
root.protocol = configure.protocol;
struct2xml(root, [rootpath 'system\mod\sample_protocol.xml']);

% clean configure
configure = configureclean(configure);

% output sample xmls 2
root = [];
root.configure = configure;
struct2xml(root, [rootpath 'system\mod\sample_output_configure.xml']);

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
    
    
    
end


