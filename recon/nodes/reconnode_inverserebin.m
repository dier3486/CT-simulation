function [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status)
% cali node, inverse of the rebin
% [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status);
% use to find out the inverse of the rebin

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

% parameters to use
Nshot = prmflow.rebin.Nshot;
focalspot = prmflow.rebin.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
Nfocal = prmflow.rebin.Nfocal;
Nviewprot_inv = prmflow.rebin.Nviewprot;
Nviewprot = Nviewprot_inv / Nfocal;
delta_view_inv = prmflow.rebin.delta_view;
delta_view = delta_view_inv * Nfocal;
% detector
detector = prmflow.system.detector;
Npixel_orig = double(detector.Npixel);
Nslice = double(detector.Nmergedslice);
% mid_U = single(detector.mid_U);
delta_d = prmflow.rebin.delta_d;  
if isfield(prmflow.recon, 'Npixel')
    Npixel_reb = prmflow.recon.Npixel;              % (could) replaced by reconnode_watergoback
else
    Npixel_reb = prmflow.rebin.Nreb;
end
if isfield(prmflow.recon, 'midchannel')
    midchannel = prmflow.recon.midchannel;          % (could) replaced by reconnode_watergoback
else
    midchannel = prmflow.rebin.midU_phi;
end
% I know that is the midchannel after rebin
startvindex = prmflow.rebin.startvindex;        % returned from reconnode_rebinprepare
% if isfield(prmflow.rebin, 'Nleft')    % deleted
%     Nleft = prmflow.rebin.Nleft;
% else
%     Nleft = 0;
% end

% NO QDO!

% fan angles & focal angle(s)
[fanangles, focalangle] = detpos2fanangles(detector.position, focalposition);

% % fan angles (copied from rebinprepare.m)
% y = detector.position(1:Npixel_orig, 2) - focalposition(:, 2)';
% x = detector.position(1:Npixel_orig, 1) - focalposition(:, 1)';
% fanangles = atan2(y, x);
% % focal angle(s)
% focalangle = atan2(-focalposition(:, 2), -focalposition(:, 1))';

% inverse rebin samples on theta-d space
fv = (fanangles - pi/2)./delta_view + (0:Nfocal-1)./Nfocal;
fv = reshape(mod(fv(:)+(0:Nviewprot-1), Nviewprot)+1, Npixel_orig, Nslice, Nviewprot_inv);
dv = detector.SID.*sin(fanangles - focalangle)./delta_d + midchannel;
dv = reshape(repmat(dv, 1, Nviewprot), Npixel_orig, Nslice, Nviewprot_inv);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel_reb, Nslice, Nviewprot, Nshot);

% get the index range
if isfield(dataflow.rawhead, 'index_range')
    index_range = reshape(dataflow.rawhead.index_range, 2, Nslice, Nviewprot, Nshot);
else
    index_range = ones(2, Nslice, Nviewprot, Nshot);
    index_range(2,:,:,:) = Npixel_reb;
end

% ini
D = zeros(Npixel_orig, Nslice, Nviewprot_inv, Nshot);
Idx = zeros(2, Nslice, Nviewprot_inv, Nshot);
% loop shots to inversew rebin the rawdata and index_range
for ishot = 1:Nshot
    % loop slices to do interp2 (low performance method)
    for islice = 1:Nslice
        % interp rawdata
        A = squeeze(dataflow.rawdata(:, islice, :, ishot));
        A = [A(:, end-startvindex+2:end) A(:, 1:end-startvindex+1)];
        A = [A A(:,1)];
        D(:, islice, :, ishot) = interp2(A, squeeze(fv(:, islice, :)), squeeze(dv(:, islice, :)), 'linear', 0);
        % interp index range
        B = zeros(Npixel_reb, Nviewprot);
        index_ii = squeeze(index_range(:, islice, :, ishot));
        for iview = 1:Nviewprot
            B(index_ii(1, iview):index_ii(2, iview), iview) = 1;
        end
        B = [B(:, end-startvindex+2:end) B(:, 1:end-startvindex+1)];
        B = [B B(:,1)];
        B = interp2(B, squeeze(fv(:, islice, :)), squeeze(dv(:, islice, :)), 'linear', 0);
        for iview = 1:Nviewprot_inv
            s = B(:, iview)>0.999;
            Idx(1, islice, iview, ishot) = find(s, 1, 'first');
            Idx(2, islice, iview, ishot) = find(s, 1, 'last');
        end
    end
end
% reshape
D = reshape(D, Npixel_orig*Nslice, []);
Idx = reshape(Idx, 2*Nslice, []);

% inverse the view angle and startviewangle
viewangle = reshape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
viewangle_inv = [viewangle(end-startvindex+2 :end, :); viewangle(1 : end-startvindex+1, :)];
% I know for multi-shots the startvindex is same for each shot, and the startviewangle are not always equal.
% startviewangle = viewangle_inv(1, :);
% for DFS
viewangle_inv = repelem(viewangle_inv, Nfocal, 1) + repmat((0:Nfocal-1)'.*delta_view_inv, Nviewprot, 1);
viewangle_inv = mod(viewangle_inv, pi*2);

% to return
dataflow.rawdata = D;
dataflow.rawhead.index_range = Idx;
dataflow.rawhead.viewangle = viewangle_inv(:)';
% NOTE: the dataflow.rawhead.refblock can not be inversed, user may call the databackup node before rebin to keep it.

% to recover the values in .recon which replaced by reconnode_watergoback
prmflow.recon.Npixel = prmflow.rebin.Nreb;
prmflow.recon.midchannel = prmflow.rebin.midU_phi;

% prmflow.recon.Npixel = Npixel_orig;
% prmflow.recon.Nviewprot = Nviewprot_inv;
% prmflow.recon.Nview = prmflow.recon.Nview*Nfocal;
% prmflow.recon.delta_view = delta_view_inv;
% prmflow.recon.startviewangle = startviewangle;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end