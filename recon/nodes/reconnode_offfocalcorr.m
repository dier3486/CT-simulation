function [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status)
% recon node, off-focal correction
% [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status);
% only for axial

% parameters of this node in pipe
if ~isempty(status)
    caliprm = prmflow.pipe.(status.nodename);
else
    % for debug
    caliprm = struct();     
end
% parameters to use in prmflow
Nview = prmflow.recon.Nview;
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nviewprot = prmflow.recon.Nviewprot;
focalspot = prmflow.system.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% detector
detector = prmflow.system.detector;

% calibration table
% if isfield(prmflow.corrtable, status.nodename)
%     offcorr = prmflow.corrtable.(status.nodename);
% else
%     offcorr = struct();
% end

% test
offcorr = prmflow.corrtable.Beamharden;

% merge the caliprm to offcorr
offcorr = structmerge(caliprm, offcorr, 1, 0);
airrate = 2.^(-offcorr.airrate);
airrate = max(airrate)./airrate;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% air rate
% dataflow.rawdata = dataflow.rawdata+offcorr.airrate;
% exp
dataflow.rawdata = 2.^(-dataflow.rawdata);
% mean on slices
Aoff = squeeze(mean(reshape(dataflow.rawdata, Npixel, Nslice, Nview), 2));
% offset fix
Aoff = offfocalconv(Aoff, detector, focalposition, Nviewprot, offcorr.offwidth, offcorr.offintensity, offcorr.offedge);
% fix rawdata (with airrate)
dataflow.rawdata = dataflow.rawdata.*(1+offcorr.offintensity.*airrate) - repmat(Aoff, Nslice, 1).*airrate;
% % orig no airrate
% dataflow.rawdata = dataflow.rawdata.*(1+offcorr.offintensity) - repmat(Aoff, Nslice, 1);

% log2
dataflow.rawdata = -log2(dataflow.rawdata);


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end