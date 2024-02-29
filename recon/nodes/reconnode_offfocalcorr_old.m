function [dataflow, prmflow, status] = reconnode_offfocalcorr_old(dataflow, prmflow, status)
% recon node, off-focal correction
% [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status);
% only for axial

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% parameters of this node in pipe
if ~isempty(status)
    caliprm = prmflow.pipe.(status.nodename);
else
    % for debug
    caliprm = struct();     
end
% parameters to use in prmflow
Nshot = prmflow.raw.Nshot;
Nview = prmflow.raw.Nview;
Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
Nviewprot = prmflow.raw.Nviewprot;
scantype = prmflow.raw.scan;
% focalspot = prmflow.recon.focalspot;
% focalposition = prmflow.system.focalposition(focalspot, :);
% fix to QFS
focalposition = prmflow.system.focalposition(1, :);
% detector
detector = prmflow.system.detector;
SID = detector.SID;
SDD = detector.SDD;

% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
% mean
fanangles = mean(reshape(fanangles, Npixel, Nslice), 2);

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

% exp
dataflow.rawdata = 2.^(-dataflow.rawdata);

% z-cross
if isfield(offcorr, 'crossrate')
    if isfinite(offcorr.crossrate)
        crossrate = offcorr.crossrate;
    else
        warning('Illeagal value of ''crossrate'' in off-focal calibration table! Replaced by 0.');
        crossrate = 0;
    end
else
    crossrate = 0;
    % 0 is the mean of all the slices
end

switch lower(scantype)
    case {'axial', 'static'}
        % Axial or static
        % reshape
        dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nviewprot, Nshot);
        for ishot = 1:Nshot
            % prepare the off-focal base intensity with a z-crossed method
            Aoff = offfocalzcross(reshape(dataflow.rawdata(:,:, ishot), Npixel, Nslice, Nviewprot), crossrate);
            % offfocalfix
            offfocalfix = offfocalconvAxial(Aoff, fanangles, SID, SDD, Nviewprot, offcorr.offwidth, ...
                offcorr.offintensity, offcorr.offedge);
            offfocalfix = reshape(offfocalfix, Npixel*Nslice, []).*airrate;
            % fix rawdata (with airrate)
            dataflow.rawdata(:,:, ishot) = dataflow.rawdata(:,:, ishot).*(1+airrate.*sum(offcorr.offintensity)) - offfocalfix;
        end
    case 'helical'
        % Helical
        % reshape
        dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
        % prepare the off-focal base intensity with a z-crossed method
        Aoff = offfocalzcross(reshape(dataflow.rawdata, Npixel, Nslice, Nview), crossrate);
        % offfocalfix
        offfocalfix = offfocalconvHelical(Aoff, fanangles, SID, SDD, Nviewprot, offcorr.offwidth, ...
            offcorr.offintensity, offcorr.offedge);
        offfocalfix = reshape(offfocalfix, Npixel*Nslice, []).*airrate;
        % fix rawdata (with airrate)
        dataflow.rawdata = dataflow.rawdata.*(1+airrate.*sum(offcorr.offintensity)) - offfocalfix;
    otherwise
        % what?
        warning('Illeagal scan type for off-focal correction: %s!', scantype);
        return;
end

% min cut
minval = 2^-32;
dataflow.rawdata(dataflow.rawdata<minval) = minval;
% log2
dataflow.rawdata = -log2(dataflow.rawdata);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end