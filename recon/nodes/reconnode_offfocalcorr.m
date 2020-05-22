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
scantype = prmflow.recon.scan;
% focalspot = prmflow.recon.focalspot;
% focalposition = prmflow.system.focalposition(focalspot, :);
% fix to QFS
focalposition = prmflow.system.focalposition(1, :);
% detector
detector = prmflow.system.detector;

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    offcorr = prmflow.corrtable.(status.nodename);
else
    error('Did not load corrtable for %s!', status.nodename);
end

% off-focal kernel parameter (used in cali)
if isfield(caliprm, 'offfocalkernel')
    offfocalkernel = caliprm.offfocalkernel;
elseif isfield(prmflow.system, 'offfocalkernel')
    offfocalkernel = prmflow.system.offfocalkernel;
else
    offfocalkernel = [];
end
% load offfocalkernel file
if ischar(offfocalkernel)
    offfocalkernel = readcfgfile(offfocalkernel);
end
offcorr_fromkernel = offfocalloadkernel(offfocalkernel, prmflow.protocol);

% merge the offcorr_fromkernel and/or caliprm to offcorr
caliprm = structmerge(caliprm, offcorr_fromkernel, 1, 0); 
offcorr = structmerge(caliprm, offcorr, 1, 0);

% multi off-focal kernel?
% Noff = size(offcorr.offintensity(:), 1);

% airrate
% airrate = repmat(2.^(-offcorr.airrate(:)), 1, Noff);
% if isfield(offcorr, 'ratescale')
%     for ii = 1:Noff
%         airrate(:, ii) = airrate(:, ii).^offcorr.ratescale(ii);
%     end
% end
airrate = 2.^(-offcorr.airrate(:).*offcorr.ratescale(1));
airrate = max(airrate, 1)./airrate;
% now we only employ one ratescale to apply on all steps 

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% air rate
% dataflow.rawdata = dataflow.rawdata+offcorr.airrate;
% exp
dataflow.rawdata = 2.^(-dataflow.rawdata);

% z-cross
if isfield(offcorr, 'crossrate')
    crossrate = offcorr.crossrate;
else
    crossrate = 0;
    % 0 is the mean of all the slices
end
% prepare the off-focal base intensity with a z-crossed method
Aoff = offfocalzcross(reshape(dataflow.rawdata, Npixel, Nslice, Nview), crossrate);

switch lower(scantype)
    case {'axial', 'static'}
        % Axial or static
        offfocalfix = offfocalconvAxial(Aoff, detector, focalposition, Nviewprot, offcorr.offwidth, ...
            offcorr.offintensity, offcorr.offedge);
    case 'helical'
        % Helical
        offfocalfix = offfocalconvHelical(Aoff, detector, focalposition, Nviewprot, offcorr.offwidth, ...
            offcorr.offintensity, offcorr.offedge);
    otherwise
        % what?
        warning('Illeagal scan type for off-focal correction: %s!', scantype);
        return;
end
offfocalfix = reshape(offfocalfix, Npixel*Nslice, []).*airrate;

% offfocalfix = offfocalfix + repmat(offfocalfix_ii, Nslice, 1).*airrate(:, ii);

% Aoff = offfocalconv(Aoff, detector, focalposition, Nviewprot, offcorr.offwidth, offcorr.offintensity, offcorr.offedge);
% fix rawdata (with airrate)
dataflow.rawdata = dataflow.rawdata.*(1+airrate.*sum(offcorr.offintensity)) - offfocalfix;
% % orig no airrate
% dataflow.rawdata = dataflow.rawdata.*(1+offcorr.offintensity) - repmat(Aoff, Nslice, 1);

% log2
dataflow.rawdata = -log2(dataflow.rawdata);


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end