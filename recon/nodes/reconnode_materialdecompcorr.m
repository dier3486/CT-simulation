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

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end

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
if ~slicezebra
    error('Only slice zebra mode now!');
end

% calibration table
mdcorr = prmflow.corrtable.(status.nodename);

% shit
mdcorr = everything2single(mdcorr);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);

if GPUonoff
    dataflow.rawdata = gpuArray(dataflow.rawdata);
    [mdcorr, HU] = ...
        putinGPU(mdcorr, HU);
end

% prepare table Z
tableZr = (mdcorr.Tablelk2r.*(mdcorr.ZB^3 - mdcorr.ZA^3)+mdcorr.ZA^3).^(1/3);

lambda = mdcorr.LHmixlambda;

Plambda = dataflow.rawdata./HU;

% LaplacedTV smooth
P0 = zeros(Npixel, Nslice+2, Nview, 'like', dataflow.rawdata);
P0(:, 2:end-1, :) = Plambda;
P0(:, 1, :) = P0(:, 3, :).*1.5 - P0(:, 5, :).*0.5;
P0(:, end, :) = P0(:, end-2, :).*1.5 - P0(:, end-4, :).*0.5;
P0 = LaplacedTV1D2(P0, 0.5, 0.1, 0.2, [], 2, [], 50, 0.01);

% kappa pair based on smoothed tendency interpolation
Pkappa = (Plambda - P0(:,2:end-1,:))./2;
Pkappa = P0(:,2:end-1,:) + Pkappa(:, [2 1:end-1], :) + Pkappa(:, [2:end end-1], :);
Pkappa = Plambda - Pkappa;
Pkappa(:, 2:2:end, :) = -Pkappa(:, 2:2:end, :);

% sign
Plsign = sign(Plambda);
Plk_sign = Plsign.*sign(Pkappa);
% index of sign
s_pos = Plk_sign > 0;
s_neg = Plk_sign < 0;
s_zero = Plk_sign == 0;

sigma = 1e-5;
% Pedge 
Pedge = tablerenorm(Plambda, mdcorr.Plrange, mdcorr.NPlsamp);
% Low
sL_neg = s_neg;
sL_neg(:, 2:2:end, :) = 0;
Pedge(sL_neg) = interp1(mdcorr.PLedge(:, 1),  Pedge(sL_neg), 'linear', 'extrap');
sL_pos = s_pos;
sL_pos(:, 2:2:end, :) = 0;
Pedge(sL_pos) = interp1(mdcorr.PLedge(:, 2),  Pedge(sL_pos), 'linear', 'extrap');
% High
sH_neg = s_neg;
sH_neg(:, 1:2:end, :) = 0;
Pedge(sH_neg) = interp1(mdcorr.PHedge(:, 1),  Pedge(sH_neg), 'linear', 'extrap');
sH_pos = s_pos;
sH_pos(:, 1:2:end, :) = 0;
Pedge(sH_pos) = interp1(mdcorr.PHedge(:, 2),  Pedge(sH_pos), 'linear', 'extrap');
% 0
Pedge(s_zero) = 0;
% P_kappa
Pkappa = hardentanh2(Pkappa, sigma, Pedge);
% P_lambda
Plambda(:, 1:2:end, :) = Plambda(:, 1:2:end, :) - Pkappa(:, 1:2:end, :).*lambda;
Plambda(:, 2:2:end, :) = Plambda(:, 2:2:end, :) + Pkappa(:, 2:2:end, :).*(1-lambda);

% post mu
mu_ref = 0.025;
ErfSigma = 0.05;
muErf = erf(sqrt(exp(-Plambda.*mu_ref)).*abs(Plambda).*ErfSigma);
muErf(muErf<0.2) = 0.2;

% renorm Pl
Plambda = tablerenorm(Plambda, mdcorr.Plrange, mdcorr.NPlsamp);

% look up Pkrange
[NPksamp_neg, NPksamp_pos] = tac(mdcorr.NPksamp);
% NPksamp = NPksamp_neg + NPksamp_pos - 1;
Pkrange = zeros(Npixel, Nslice, Nview, 'like', dataflow.rawdata);
Pkrange(s_neg) = 1./interp1(mdcorr.Pkrange(:, 1), Plambda(s_neg));
Pkrange(s_pos) = 1./interp1(mdcorr.Pkrange(:, 2), Plambda(s_pos));

% renorm Pk
Pkappa = Pkappa.*Plsign.*Pkrange;
Pkappa = Pkappa.*(s_neg.*(1 - NPksamp_neg) + s_pos.*(NPksamp_pos - 1)) + NPksamp_neg;

% % boundary (shall need not)
% Pkappa(Pkappa<1) = 1;
% Pkappa(Pkappa>NPksamp) = NPksamp;

% look up table Z
R = interp2(tableZr, Pkappa, Plambda);
% denoise
for islice = 1:Nslice
    R_ii = squeeze(R(:, islice, :));
    mu_ii = squeeze(muErf(:, islice, :));
    R(:, islice, :) = BregmanTV2D(R_ii, mu_ii, mu_ii./5);
end

% R
R = (R.^3 - mdcorr.ZA^3)./(mdcorr.ZB^3 - mdcorr.ZA^3);

% renorm R
[NRsamp_neg, NRsamp_pos] = tac(mdcorr.NRsamp);
% NRsamp = NRsamp_neg + NRsamp_pos - 1;
s_neg = R < 0;
Rnorm = R.*(s_neg.*((1 - NRsamp_neg)/mdcorr.Rrange(1)) + ~s_neg.*((NRsamp_pos - 1)/mdcorr.Rrange(2))) + NRsamp_neg;

% look up D
D = interp2(mdcorr.Tablelr2D, Rnorm, Plambda);
D = D.*Plsign;

% return rawdata
dataflow.rawdata = gather(reshape(D.*(1+R.*1i).*HU, Npixel*Nslice, Nview));

end

function Pout = tablerenorm(Pin, range, Nsamp)

Pout = (Pin - range(1)) ./ (range(2) - range(1));
Pout = abs(Pout).*(Nsamp - 1) + 1;
Pout(Pout > Nsamp) = Nsamp;

end


    