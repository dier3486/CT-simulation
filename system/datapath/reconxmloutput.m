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
    % only Axial supported yet
    recon{iw}.pipe.Axialrebin = struct();
    % QDO rebin
    recon{iw}.pipe.Axialrebin.QDO = 0;
    % hard code FBP for temprory use
%     recon{iw}.pipe.FBP = struct();
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

% detector corr table
system.detector_corr = SYS.detector.detector_corr.frame_base;
% focal position(s)
system.focalposition = SYS.source.focalposition;
% DCB
if isfield(SYS, 'datacollector')
    system.angulationcode = SYS.datacollector.angulationcode;
    system.angulationzero = SYS.datacollector.angulationzero;
    system.DBBzero = SYS.datacollector.DBBzero;
end
% console
if isfield(SYS, 'console')
    % how the console explain the protocol
    if isfield(SYS.console.protocaltrans, 'collimatorexplain_file')
        system.collimatorexplain = SYS.console.protocaltrans.collimatorexplain_file;
    elseif isfield(SYS.console.protocaltrans, 'collimatorexplain')
        system.collimatorexplain = SYS.console.protocaltrans.collimatorexplain;
    end
    % filename rule
    if isfield(SYS.console.protocaltrans, 'filetagsrule_file')
        system.filetagsrule = SYS.console.protocaltrans.filetagsrule_file;
    elseif isfield(SYS.console.protocaltrans, 'filetagsrule')
        system.filetagsrule = SYS.console.protocaltrans.filetagsrule;
    end
    % corr couple rule
    if isfield(SYS.console.protocaltrans, 'corrcouplerule')
        system.corrcouplerule = SYS.console.protocaltrans.corrcouplerule;
    end
    % nominal slice thickness
    if isfield(SYS.console.protocaltrans, 'nominalslicethickness')
        system.nominalslicethickness = SYS.console.protocaltrans.nominalslicethickness;
    end
end
end