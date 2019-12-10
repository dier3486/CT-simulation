function simuresultsoutput(SYS, Data)
% simuresultsoutput(SYS, Data)

% output the rawdata and air (no offset?)
rawdataoutput(SYS, Data);

% to be use
% mu_ref = 0.020323982562342;
mu_ref = 0.021139124532511;
HCscale = 1000*log(2)/mu_ref;

system = systemforrecon(SYS);

% recon xml
Nw = SYS.source.Wnumber;
recon = cell(1, Nw);
for iw = 1:Nw
    % rawdata
    recon{iw}.rawdata = [SYS.output.path SYS.output.files.rawdata{iw} '_' SYS.output.rawdataversion '.raw'];
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
    recon{iw}.pipe.Air.corr = [SYS.output.path SYS.output.files.aircorr{iw} '_v1.0' '.corr'];
    recon{iw}.pipe.Beamharden = struct();
    recon{iw}.pipe.Housefield = struct();
    recon{iw}.pipe.Housefield.HCscale = HCscale;
    recon{iw}.pipe.Rebin = struct();
    recon{iw}.pipe.Filter = struct();
    recon{iw}.pipe.Backprojection = struct();
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
end