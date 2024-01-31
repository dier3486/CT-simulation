function [dataflow, prmflow, status] = reconnode_materialdecompcorr(dataflow, prmflow, status)
% recon node, two-material decomposition by low-high energy
% [dataflow, prmflow, status] = reconnode_materialdecompcorr(dataflow, prmflow, status);

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
Nview = prmflow.raw.Nview;
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
HU = 1000;
if isfield(prmflow.system, 'slicezebra')
    slicezebra = prmflow.system.slicezebra;
else
    slicezebra = false;
end

% calibration table
mdcorr = prmflow.corrtable.(status.nodename);

% shit
mdcorr = everything2single(mdcorr);

lambda = mdcorr.LHmixlambda;
Mlambda = [1-lambda lambda; 1 -1];
Minvlam = [1 lambda; 1 lambda-1];

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

% slice zebra
if slicezebra
    % to interp1
    Nslice2 = Nslice;
    % linear interp
    Plk = zeros(Npixel, Nslice, Nview, 2, 'like', dataflow.rawdata);
    % low
    Plk(:, 1:2:end, :, 1) = dataflow.rawdata(:, 1:2:end, :);
    Plk(:, 2:2:end-2, :, 1) = (dataflow.rawdata(:, 1:2:end-2, :) + dataflow.rawdata(:, 3:2:end, :))./2;
    Plk(:, end, :, 1) = dataflow.rawdata(:, end-1, :);
    % high
    Plk(:, 2:2:end, :, 2) = dataflow.rawdata(:, 2:2:end, :);
    Plk(:, 3:2:end, :, 2) = (dataflow.rawdata(:, 2:2:end-2, :) + dataflow.rawdata(:, 4:2:end, :))./2;
    Plk(:, 1, :, 2) = dataflow.rawdata(:, 2, :);
    % reshape
    Plk = reshape(Plk, [], 2);
    % 1/1000
    Plk = Plk./HU;
    % edge
    % TBC
else
    Nslice2 = Nslice/2;
    Plk = cat(2, reshape(dataflow.rawdata(:, 1:2:end, :), [], 1), reshape(dataflow.rawdata(:, 2:2:end, :), [], 1));
    % 1/1000
    Plk = Plk./HU;
end
Np = size(Plk, 1);

% L/H to lamda/kappa
Plk = Plk * Mlambda';

% renorm Pl
Plk_sign = sign(Plk);
Plk(:, 1) = (Plk(:, 1) - mdcorr.Plrange(1)) ./ (mdcorr.Plrange(2) - mdcorr.Plrange(1));
Plk(:, 1) = abs(Plk(:, 1)).*(mdcorr.NPlsamp - 1) + 1;
Plk(Plk(:, 1) > mdcorr.NPlsamp, 1) = mdcorr.NPlsamp;

% look up Pkrange
[NPksamp_neg, NPksamp_pos] = tac(mdcorr.NPksamp);
NPksamp = NPksamp_neg + NPksamp_pos - 1;
Pkrange = zeros(Np, 2, 'like', Plk);
Pkrange(:, 1) = interp1(mdcorr.Pkrange(:, 1), Plk(:, 1));
Pkrange(:, 2) = interp1(mdcorr.Pkrange(:, 2), Plk(:, 1));
% renorm Pk
s_pos = Plk_sign(:, 1) .* Plk_sign(:, 2) > 0;
Plk(s_pos, 2) = Plk(s_pos, 2)./Pkrange(s_pos, 2).*(NPksamp_pos - 1) + NPksamp_neg;
s_neg = Plk_sign(:, 1) .* Plk_sign(:, 2) < 0;
Plk(s_neg, 2) = Plk(s_neg, 2)./Pkrange(s_neg, 1).*(1 - NPksamp_neg) + NPksamp_neg;
s_zero = Plk_sign(:, 1) .* Plk_sign(:, 2) == 0;
Plk(s_zero, 2) = NPksamp_neg;
% boundary
Plk(Plk(:,2)<1, 2) = 1;
Plk(Plk(:,2)>NPksamp, 2) = NPksamp;

% look up R
R = interp2(mdcorr.Tablelk2r, Plk(:,2), Plk(:,1));

% R denoise
% TBC

% renorm R
[NRsamp_neg, NRsamp_pos] = tac(mdcorr.NRsamp);
NRsamp = NRsamp_neg + NRsamp_pos - 1;
s_pos = R >= 0;
Rnorm = R.*0;
Rnorm(s_pos) = R(s_pos)./mdcorr.Rrange(2).*(NRsamp_pos - 1) + NRsamp_neg;
s_neg = ~s_pos;
Rnorm(s_neg) = R(s_neg)./mdcorr.Rrange(1).*(-NRsamp_neg + 1) + NRsamp_neg;
% boundary
Rnorm(Rnorm<1) = 1;
Rnorm(Rnorm>NRsamp) = NRsamp;

% look up D
D = interp2(mdcorr.Tablelr2D, Rnorm, Plk(:,1));
D = D.*Plk_sign(:, 1);

% return rawdata
dataflow.rawdata = reshape(D.*(1+R.*1i), Npixel*Nslice2, Nivew);

% % test
% Plk = reshape(Plk, Npixel, Nslice, Nview, 2);
% Pkrange = reshape(Pkrange, Npixel, Nslice, Nview, 2);
% R = reshape(R, Npixel, Nslice, Nview);
% D = reshape(D, Npixel, Nslice, Nview);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end