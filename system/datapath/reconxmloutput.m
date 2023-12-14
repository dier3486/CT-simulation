function [recon, reconxmlfile] = reconxmloutput(SYS, tofile)
% output recon xml, return and/or to file
% [recon, reconxml] = reconxmloutput(SYS);
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
        warn('Unknown style %s to save the raw data!', SYS.output.rawdatastyle);
        rawext = '';
end

% make 
for iw = 1:Nw
    % rawdata
    recon{iw}.rawdata = fullfile(SYS.output.path, [SYS.output.files.rawdata{iw} rawext]);
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
    recon{iw}.pipe.Log2 = struct();
    recon{iw}.pipe.Air = struct();
    if isfield(SYS.output.files, 'air')
        recon{iw}.pipe.Air.corr = fullfile(SYS.output.path, [SYS.output.files.air{iw} '.corr']);
    end
    recon{iw}.pipe.Beamharden = struct();
    if isfield(SYS.output.files, 'beamharden')
        recon{iw}.pipe.Beamharden.corr = fullfile(SYS.output.path, [SYS.output.files.beamharden{iw} '.corr']);
    end
    recon{iw}.pipe.Hounsefield = struct();
    recon{iw}.pipe.Hounsefield.HCscale = HCscale;
    switch lower(SYS.protocol.scan)
        case 'axial'
            % Axial
            recon{iw}.pipe.Axialrebin = struct();
            % no QDO
            recon{iw}.pipe.Axialrebin.QDO = 0;
        case 'helical'
            % Helical
            recon{iw}.pipe.Helicalrebin = struct();
    end
    % filter
    recon{iw}.pipe.Filter = struct();
    recon{iw}.pipe.Filter.name = 'hann';
    recon{iw}.pipe.Filter.freqscale = 1.2;
    % BP
    recon{iw}.pipe.Backprojection = struct();
    recon{iw}.pipe.Backprojection.FOV = 500;

    % TBC
end
% save xml file
if tofile
    root.configure.recon = recon;
    reconxmlfile = fullfile(SYS.output.path, [SYS.output.files.reconxml '.xml']);
    struct2xml(root, reconxmlfile);
else
    reconxmlfile = [];
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
% DCB
if isfield(SYS, 'datacollector')
    system.angulationcode = SYS.datacollector.angulationcode;
    system.angulationzero = SYS.datacollector.angulationzero;
    system.DBBzero = SYS.datacollector.DBBzero;
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
end