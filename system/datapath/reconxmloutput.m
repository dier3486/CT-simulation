function [recon, reconcfgfile] = reconxmloutput(SYS, tofile)
% output recon .xml (or .json), return and/or to file
% [recon, reconcfgfile] = reconxmloutput(SYS);
% or recon = reconxmloutput(SYS, 0); to avoid writing to file.

if nargin<2
    tofile = true;
end

% recon.system
system = systemforrecon(SYS);

% I know
HCscale = 1000;

% ini
Nw = SYS.source.Wnumber;
recon = cell(1, Nw);
recon(:) = {struct()};

% output type
switch SYS.output.rawdatastyle
	case {'24bit', '16bit', 'single'}
        rawext = '.raw';
    case 'mat'
        rawext = '.mat';
    otherwise
        warning('Unknown style %s to save the raw data!', SYS.output.rawdatastyle);
        rawext = '';
end

% pipelinereplicate
if isfield(SYS.protocol, 'pipelinereplicate')
    pipelinereplicate = SYS.protocol.pipelinereplicate;
else
    pipelinereplicate = false;
end

% make 
for iw = 1:Nw
    % rawdata
    recon{iw}.rawdata = pathclbs(fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} rawext]));
    % IOpath
    recon{iw}.IOstandard = SYS.path.IOstandard;
    % system
    recon{iw}.system = system;
    % protocol
    recon{iw}.protocol = SYS.protocol;
    recon{iw}.protocol.KV = SYS.source.KV{iw};
    recon{iw}.protocol.mA = SYS.source.mA{iw};
    % protocol.focalspot
    recon{iw}.protocol.focalspot = sum(2.^(SYS.protocol.focalspot-1));
    % recon work flow
    % corrections
    recon{iw}.pipe.Log2 = struct();
    recon{iw}.pipe.Air = struct();
    if isfield(SYS.output.files, 'air')
        recon{iw}.pipe.Air.corr = pathclbs(fullfile(SYS.output.path, [SYS.output.files.air{iw} '.corr']));
    end
    recon{iw}.pipe.Beamharden = struct();
    if isfield(SYS.output.files, 'beamharden')
        recon{iw}.pipe.Beamharden.corr = pathclbs(fullfile(SYS.output.path, [SYS.output.files.beamharden{iw} '.corr']));
    end
    recon{iw}.pipe.Hounsfield = struct();
    recon{iw}.pipe.Hounsfield.HCscale = HCscale;

    % Rebin
    recon{iw}.pipe.Rebin = struct();

    % filter
    if ~pipelinereplicate
        recon{iw}.pipe.Filter = struct();
        recon{iw}.pipe.Filter.name = 'hann';
        recon{iw}.pipe.Filter.freqscale = 1.2;
    end
    % BP
    recon{iw}.pipe.Backprojection = struct();
    if pipelinereplicate
        % put filter in BP
        recon{iw}.pipe.Backprojection.Filter = struct();
        recon{iw}.pipe.Backprojection.Filter.name = 'hann';
        recon{iw}.pipe.Backprojection.Filter.freqscale = 1.2;
    end
    if isfield(SYS.protocol, 'reconFOV') && ~isempty(SYS.protocol.reconFOV)
        recon{iw}.pipe.Backprojection.FOV = SYS.protocol.reconFOV;
    else
        recon{iw}.pipe.Backprojection.FOV = 500;
    end
    if isfield(SYS.detector, 'concyclic') && ~isempty(SYS.detector.concyclic)
        recon{iw}.pipe.Backprojection.concyclic = SYS.detector.concyclic;
    end
    % output or pipelinestuck
    if pipelinereplicate
        recon{iw}.pipe.pipelinestuck = struct();
    end

    % TBC
end
% save .xml or .json file
if tofile
    if isfield(SYS.console, 'reconcfgfile')
        cfgext = SYS.console.reconcfgfile;
    else
        cfgext = '.xml';
    end
    root.configure.recon = recon;
    reconcfgfile = fullfile(SYS.output.path, [SYS.output.files.reconxml cfgext]);
    switch cfgext
        case '.xml'
            struct2xml(root, reconcfgfile);
        case '.json'
            jsonwrite(root, reconcfgfile);
        case '.mat'
            save(reconcfgfile, '-struct', 'root');
        otherwise
            % ??
            0;
    end
else
    reconcfgfile = [];
end

end


function system = systemforrecon(SYS)
% system paramter and data for recon

if isfield(SYS, 'console') && isfield(SYS.console, 'reconsystem')
    % congifured recon system
    system = SYS.console.reconsystem;
    return;
end

% detector corr table
system.detector_corr = SYS.detector.detector_corr.frame_base;
% focal position(s)
if isfield(SYS.source, 'tube_corr')
    system.focalposition = SYS.source.tube_corr.focalposition;
else
    system.focalposition = SYS.source.focalposition;
end
% maxFOV
if isfield(SYS.detector, 'maxFOV')
    system.maxFOV = SYS.detector.maxFOV;
else
    system.maxFOV = 500;
end
% slicezebra and ZebraOrder
if isfield(SYS.detector, 'slicezebra')
    system.slicezebra = SYS.detector.slicezebra;
    if isfield(SYS.detector, 'ZebraOrder')
        system.ZebraOrder = SYS.detector.ZebraOrder;
    else
        system.ZebraOrder = 1;
    end
else
    system.slicezebra = false;
end
% concyclic
if isfield(SYS.detector, 'concyclic')
    system.concyclic = SYS.detector.concyclic;
else
    system.concyclic = false;
end

% DCB
if isfield(SYS, 'datacollector')
    system.angulationcode = SYS.datacollector.angulationcode;
    system.angulationzero = SYS.datacollector.angulationzero;
    system.DBBzero = SYS.datacollector.DBBzero;
    system.movercode = SYS.datacollector.movercode;
    system.moverlength = SYS.datacollector.moverlength;
end
% console
if isfield(SYS, 'console')
    % how the console explain the protocol
    if isfield(SYS.console.protocoltrans, 'collimatorexplain_file')
        system.collimatorexplain = SYS.console.protocoltrans.collimatorexplain_file;
    elseif isfield(SYS.console.protocoltrans, 'collimatorexplain')
        system.collimatorexplain = SYS.console.protocoltrans.collimatorexplain;
    end
    % filename rule
    if isfield(SYS.console.protocoltrans, 'filetagsrule_file')
        system.filetagsrule = SYS.console.protocoltrans.filetagsrule_file;
    elseif isfield(SYS.console.protocoltrans, 'filetagsrule')
        system.filetagsrule = SYS.console.protocoltrans.filetagsrule;
    end
    % corr couple rule
    if isfield(SYS.console.protocoltrans, 'corrcouplerule')
        system.corrcouplerule = SYS.console.protocoltrans.corrcouplerule;
    end
    % nominal slice thickness
    if isfield(SYS.console.protocoltrans, 'nominalslicethickness')
        system.nominalslicethickness = SYS.console.protocoltrans.nominalslicethickness;
    end
end

% label
if isfield(SYS, 'label')
    % Manufacturer  (0008,0070)
    if isfield(SYS.label, 'Manufacturer')
        system.Manufacturer = SYS.label.Manufacturer;
    end
    % (Manufacturer) ModelName  (0008,1090)
    if isfield(SYS.label, 'ModelName')
        system.ModelName = SYS.label.ModelName;
    end
    % DeviceSerialNumber  (0018,1000)
    if isfield(SYS.label, 'DeviceSerialNumber') 
        system.DeviceSerialNumber = SYS.label.DeviceSerialNumber;
    end
end
end