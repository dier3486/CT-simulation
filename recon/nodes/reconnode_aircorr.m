function [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);

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

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nview = prmflow.recon.Nview;
Nfocal = prmflow.recon.Nfocal;
% parameters set in pipe
airprm = prmflow.pipe.(status.nodename);

% calibration table
aircorr = prmflow.corrtable.(status.nodename);
% parameters in corr
Nsect = double(aircorr.Nsection);
refpixel = double(aircorr.refpixel);
Nref = double(aircorr.refnumber);
if isfield(aircorr, 'referrcut')
    referrcut = reshape(aircorr.referrcut, Nref, Nfocal);
else
    referrcut = zeros(Nref, Nfocal);
end
% slice independent refernece?
if isfield(airprm, 'sliceindependent')
    sliceindependent = airprm.sliceindependent;
else
    sliceindependent = false;
end

% angles of the air table
sectangle = (pi*2/Nsect);
% airangle = (-1:Nsect).*(pi*2/Nsect);
% airmain & airref
aircorr.main = reshape(aircorr.main, [], Nsect);
airmain = [aircorr.main aircorr.main(:,1)];
if isfield(aircorr, 'referenceKVmA')
    airKVmA = reshape(aircorr.referenceKVmA, [], Nsect);
    if size(airKVmA, 1) < Nfocal
        % This should be a bug in air calibration table
        airKVmA = repmat(airKVmA, Nfocal, 1);
    end
    airKVmA = [airKVmA airKVmA(:, 1)];
else
    airKVmA = zeros(Nfocal, Nsect+1);
end

% interpolation index and weight
retangle = mod(dataflow.rawhead.viewangle - aircorr.firstangle, pi*2);
intp_index = floor(retangle./sectangle);
intp_alpha = retangle./sectangle - intp_index;
intp_index = intp_index + 1;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

% KVmA
KVmA = -log2(dataflow.rawhead.KV.*dataflow.rawhead.mA);

% corr rawdata with air
for ifocal = 1:Nfocal
    % rawdata
    viewindex = ifocal:Nfocal:Nview;
    airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal-1);
    dataflow.rawdata(:, viewindex) = ...
        dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index(viewindex)).*(1-intp_alpha(viewindex));
    dataflow.rawdata(:, viewindex) = ...
        dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index(viewindex)+1).*intp_alpha(viewindex);
    % KVmA
    KVmA(viewindex) = ...
        KVmA(viewindex) - airKVmA(ifocal, intp_index(viewindex)).*(1-intp_alpha(viewindex)) ...
        - airKVmA(ifocal, intp_index(viewindex)+1).*intp_alpha(viewindex);
end

% rawref
[rawref, referr] = airreference2(dataflow.rawdata, refpixel, Npixel, Nslice, sliceindependent);

% to save refblock in rawhead
dataflow.rawhead.refblock = false(Nref, Nview);

% corr rawdata with ref
for ifocal = 1:Nfocal
    viewindex = ifocal:Nfocal:Nview;
    [rawref_ifc, refblock] = referenceblock(rawref(:, viewindex), referr(:, viewindex), ...
                 referrcut(:, ifocal), KVmA(viewindex), airprm);
    if sliceindependent
        dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - repelem(rawref_ifc, Npixel, 1);
    else
        dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - rawref_ifc;
    end
    dataflow.rawhead.refblock(:, viewindex) = refblock;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [ref, refblk] = referenceblock(rawref, referr, referrcut, KVmA, airprm)
% reference block

if isfield(airprm, 'cutscale')
    referrcut = referrcut.*airprm.cutscale;
else
    % default cut scale is 1.2;
    referrcut = referrcut.*1.2;
end

if isfield(airprm, 'blockwindow')
    m_blk = airprm.blockwindow;
else
    m_blk = 10;
end
blk1 = conv(double(referr(1,:)>max(referrcut(1), 1e-5)), ones(1, 2*m_blk+1));
blk2 = conv(double(referr(2,:)>max(referrcut(2), 1e-5)), ones(1, 2*m_blk+1));
blk1 = blk1(m_blk+1:end-m_blk)>0;
blk2 = blk2(m_blk+1:end-m_blk)>0;
idx_both = ~blk1 & ~blk2;
refblk = [blk1; blk2];

Nview = size(rawref, 2);
rawref = reshape(rawref, [], 2, Nview);
Nref = size(rawref, 1);
ref = zeros(Nref, Nview);
ref(:, idx_both) = squeeze((rawref(:, 1, idx_both) + rawref(:, 2, idx_both))./2);

if any(~idx_both)
    % some views are blocked
    blkidx_1 = ~blk1 & blk2;
    blkidx_2 = blk1 & ~blk2;
    %
    view_bkl1 = find(~idx_both, 1, 'first');
    view0 = find(idx_both, 1, 'first');
    if isempty(view0)
        % all views are blocked ??
        ref = repmat(KVmA, Nref, 1);
        return;
    end
    if view_bkl1==1
        %1 go back
        for ii = view0-1:-1:1
            if blkidx_1(ii)
                % ref1
                ref(:, ii) = ref(:, ii+1)-squeeze(rawref(:, 1, ii+1)-rawref(:, 1, ii));
            elseif blkidx_2(ii)
                % ref2
                ref(:, ii) = ref(:, ii+1)-squeeze(rawref(:, 2, ii+1)-rawref(:, 2, ii));
            else
                % mA
                ref(:, ii) = ref(:, ii+1)-KVmA(ii+1)+KVmA(ii);
            end
        end
        %2 forward
        view_bkl1 = find(~idx_both(view0:end), 1, 'first');
    end
    for ii = view_bkl1:Nview
        if idx_both(ii)
            continue
        end
        if blkidx_1(ii)
            % ref1
            ref(:, ii) = ref(:, ii-1)-squeeze(rawref(:, 1, ii-1)-rawref(:, 1, ii));
        elseif blkidx_2(ii)
            % ref2
            ref(:, ii) = ref(:, ii-1)-squeeze(rawref(:, 2, ii-1)-rawref(:, 2, ii));
        else
            % mA
            ref(:, ii) = ref(:, ii-1)-KVmA(ii-1)+KVmA(ii);
        end
    end
end

end