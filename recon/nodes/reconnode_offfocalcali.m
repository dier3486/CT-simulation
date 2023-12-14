function [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status)
% recon node, off-focal calibration based on scanning a block column 
% [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status);
%

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

% parameters set in pipe
caliprm = prmflow.pipe.(status.nodename);

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v2.0';
end

Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
Nviewprot = prmflow.raw.Nviewprot;
Nview = size(dataflow.rawdata, 2);

if isfield(caliprm, 'offsample')
    Noffsample = caliprm.offsample;
else
    Noffsample = max(2^ceil(log2(Npixel)), 64);  % =1024
end
if isfield(caliprm, 'tauremeasure')
    at3_scale = caliprm.tauremeasure.t3;
else
    at3_scale = 0.1;
end

% special data correction
% I know the air calibration table is
aircorr = prmflow.corrtable.Air;
airmean = mean(reshape(aircorr.main, [], aircorr.Nsection), 2);
minttime = mean(dataflow.rawhead.Integration_Time);

airfix = 2.^(-airmean + log2(minttime));
raw0 = 2.^(-dataflow.rawdata - airmean + log2(minttime));
airfix = reshape(airfix, Npixel, Nslice);
raw0 = reshape(raw0, Npixel, Nslice, Nview);

% L-H 
raw0 = (raw0(:, 1:2:end, :) - raw0(:, 2:2:end, :))./(airfix(:, 1:2:end)-airfix(:, 2:2:end));

% % offset corr
% dataflow.rawdata = dataflow.rawdata - mean(dataflow.offset.rawdata, 2);
% % air corr
% % I know the air calibration table is
% aircorr = prmflow.corrtable.(status.nodename);
% Nsect = single(aircorr.Nsection);
% aircorr.main = reshape(aircorr.main, [], Nsect);
% airmain = [aircorr.main aircorr.main(:,1)];
% 
% minttime = mean(dataflow.rawhead.Integration_Time);
% airmain = 2.^(-airmain+log2(minttime));
% 
% retangle = mod(dataflow.rawhead.viewangle - aircorr.firstangle, pi*2)./(pi*2/Nsect) + 1;
% airfix = interp1(airmain', retangle)';
% 
% airfix = reshape(airfix, Npixel, Nslice, Nview);
% dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
% 
% % L-H 
% raw0 = (dataflow.rawdata(:, 1:2:end, :) - dataflow.rawdata(:, 2:2:end, :))./(airfix(:, 1:2:end, :)-airfix(:, 2:2:end, :));

% focalposition
focalspot = prmflow.raw.focalspot;
focalposition = prmflow.system.focalposition(focalspot, :);
% no DFS!

% detector
detector = prmflow.system.detector;

% fanangles
[fanangles, ~] = detpos2fanangles(detector.position, focalposition);
fanangles = reshape(fanangles, Npixel, Nslice);

% half slice
Nslice2 = Nslice/2;
fanangles = (fanangles(:, 1:2:end) + fanangles(:, 2:2:end))./2;


SID = detector.SID;
SDD = detector.SDD;
alpha = acos(SID/SDD);
phi = fanangles-pi/2;
phi_off = phi - atan( sin(phi).*sin(alpha)./(cos(phi)-cos(alpha)) ) ./ sin(alpha);

t_off = double(phi - phi_off);
dt = (max(t_off(:)) - min(t_off(:)))/(Noffsample - 1);
t0 = linspace(min(t_off(:)), max(t_off(:)), Noffsample)';

% tau-measure scale
tb = max(abs(t_off(:)));
t_resp = (t0 + t0.^3.*at3_scale).*(tb/(tb + tb^3*at3_scale));
dt_resp = (1 + t0.^2.*at3_scale.*3).*(tb/(tb + tb^3*at3_scale));
dt_resp_off = (1 + t_off.^2.*at3_scale.*3).*(tb/(tb + tb^3*at3_scale));

index_p = single(1:Npixel)';
index_v = 1:Nviewprot;
delta_view = pi*2/Nviewprot;
Dphi = phi - phi_off;

% t1 = interp1(t_off, index_p, t_resp(:), 'linear', 'extrap');
raw1 = zeros(Noffsample, Nslice2, Nviewprot);
for islice = 1:Nslice2
    Df = mod(index_v - Dphi(:, islice)./delta_view - 1, Nviewprot) + 1;
    raw0_isl = [squeeze(raw0(:, islice, :)) raw0(:, islice, 1)];
    rawtmp = interp2(raw0_isl , Df, repmat(index_p, 1, Nviewprot));

    raw1(:, islice, :) = interp1(t_off(:, islice), rawtmp, t_resp(:), 'linear', 'extrap');
end

raw2 = squeeze(mean(raw1, 2));
Imin = mean(min(raw2,[],1));

data0 = ones(Noffsample, Nviewprot);
data_d = zeros(Noffsample, Nviewprot);

m = 40;
mk = 20;
u0 = 2.5;
options = optimoptions('lsqnonlin','Display','off');
for iview = 1:Nviewprot
    index_l = find(raw2(:, iview)<0.5, 1, 'first');
    index_r = find(raw2(:, iview)<0.5, 1, 'last');
    indexLR = (max(index_l-m,1) : min(index_r+m, Noffsample))';
    data_r = raw2(indexLR, iview);
    umin = mean(mink(data_r, mk));
    u = lsqnonlin(@(u) offwell(indexLR, u, umin) - data_r, [index_l, index_r, 1, 1], [], [], options);
    ufix = u;
    ufix([3 4]) = u([3 4]).*u0;
    Ru = offwell(indexLR, u, umin);
    Rfix = offwell(indexLR, ufix, umin);
    data0(indexLR, iview) = Rfix;
    data_d(indexLR, iview) = data_r - Rfix;
end

1;
% offfocalcorr
offfocalcorr = caliprmforcorr(prmflow, corrversion);
% merge
offfocalcorr = structmerge(offfocalcorr, offcorrprm);
offfocalcorr = structmerge(offfocalcorr, offfocalbase);

% steps
if isfield(offfocalcorr, 'offintensity')
    offfocalcorr.stepnumber = length(offfocalcorr.offintensity);
else
    offfocalcorr.stepnumber = 1;
end

% to return 
dataflow.offfocalcorr =  offfocalcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

function y = offwell(x, u, umin)
y = ( 1./(1+exp((x-u(1)).*u(3))) + 1./(1+exp(-(x-u(2)).*u(4))) ).*(1-umin) + umin;

end