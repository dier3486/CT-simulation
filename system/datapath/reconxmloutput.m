function recon = reconxmloutput(SYS)
% output recon xml
% recon = reconxmloutput(SYS);


system = systemforrecon(SYS);
% I know
HCscale = 1000;

% recon xml
Nw = SYS.source.Wnumber;
recon = cell(1, Nw);
recon(:) = {struct()};
for iw = 1:Nw
    % rawdata
    recon{iw}.rawdata = [SYS.output.path SYS.output.files.rawdata{iw} '.raw'];
    % IOpath
    recon{iw}.IOstandard = SYS.path.IOstandard;
    % system
    recon{iw}.system = system;
    % protocol
    recon{iw}.protocol = SYS.protocol;
    recon{iw}.protocol.KV = SYS.source.KV{iw};
    recon{iw}.protocol.mA = SYS.source.mA{iw};
    % recon work flow
    recon{iw}.pipe.Air = struct();
    if isfield(SYS.output.files, 'air')
        recon{iw}.pipe.Air.corr = [SYS.output.path SYS.output.files.air{iw} '.corr'];
    end
    recon{iw}.pipe.Beamharden = struct();
    if isfield(SYS.output.files, 'beamharden')
        recon{iw}.pipe.Beamharden.corr = [SYS.output.path SYS.output.files.beamharden{iw} '.corr'];
    end
    recon{iw}.pipe.Housefield = struct();
    recon{iw}.pipe.Housefield.HCscale = HCscale;
    % only Axial supported yet
    recon{iw}.pipe.Axialrebin = struct();
    % QDO rebin
    recon{iw}.pipe.Axialrebin.QDO = 1;
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
root.configure.recon = recon;
reconxmlfile = [SYS.output.path 'recon_series' num2str(SYS.protocol.series_index) '.xml'];
struct2xml(root, reconxmlfile);

end

function system = systemforrecon(SYS)
% system paramter and data for recon
system.detector_corr = SYS.detector.detector_corr.frame_base;
system.focalposition = SYS.source.focalposition;
if isfield(SYS, 'datacollector')
    system.angulationcode = SYS.datacollector.angulationcode;
    system.angulationzero = SYS.datacollector.angulationzero;
    system.DBBzero = SYS.datacollector.DBBzero;
end
% TBC
end