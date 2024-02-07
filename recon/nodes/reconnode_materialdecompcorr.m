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
    % 1/1000
    Plk = Plk./HU;
    
    % Pkappa denoise
    Plk_sign = sign(Plk);
    Plamda0 = (Plk(:,:,:, 1).*(1 - lambda) + Plk(:,:,:, 2).*lambda);
    Pkappa0 = (Plk(:,:,:, 1) - Plk(:,:,:, 2)) ./ (abs(Plamda0) + 2.0);
    Pkappa0 = permute(Pkappa0, [1 3 2]);

    Pkappa = BregmanTV3D_tanh(Pkappa0, 0.5, 0.1, Pkappa0, [], [], [], 2);
%     Pkappa = BregmanTV3D(Pkappa0, 0.5, 0.1, Pkappa0, [], [], [], 2);
    Pkappa = permute(Pkappa, [1 3 2]);
    Pkappa = Pkappa .* (abs(Plamda0) + 2.0);

%     Pkappa = Pkappa0.*0;
%     for islice = 1:Nslice2
%         Pkappa(:, islice, :) = BregmanTV2D(squeeze(Pkappa0(:, islice, :)), 0.5, 0.1);
%     end
    
    % edge

    % renorm PL
    PL = tablerenorm(Plk(:, 1:2:end, :, 1), mdcorr.Plrange, mdcorr.NPlsamp);
    % PLedge
    PLedge1 = interp1(mdcorr.PLedge(:, 1),  PL, 'linear', 'extrap');
    PLedge2 = interp1(mdcorr.PLedge(:, 2),  PL, 'linear', 'extrap');
    
    sigma = 1e-5;
    % Pk in PLedge
    sign_L = Plk_sign(:, 1:2:end, :, 1).*Plk_sign(:, 1:2:end, :, 2);
    PL_k = Pkappa(:, 1:2:end, :);
    PL_k(sign_L<0) = hardentanh2(PL_k(sign_L<0), sigma, PLedge1(sign_L<0));
    PL_k(sign_L>0) = hardentanh2(PL_k(sign_L>0), sigma, PLedge2(sign_L>0));
    PL_k(sign_L==0) = 0;
    Plk(:, 1:2:end, :, 1) = Plk(:, 1:2:end, :, 1) - PL_k.*lambda;
    Plk(:, 1:2:end, :, 2) = PL_k;
    
    % renorm PH
    PH = tablerenorm(Plk(:, 2:2:end, :, 2), mdcorr.Plrange, mdcorr.NPlsamp);
    % PHedge
    sign_H = Plk_sign(:, 2:2:end, :, 1).*Plk_sign(:, 2:2:end, :, 2);
    PHedge1 = interp1(mdcorr.PHedge(:, 1),  PH(sign_H<0), 'linear', 'extrap');
    PHedge2 = interp1(mdcorr.PHedge(:, 2),  PH(sign_H>0), 'linear', 'extrap');
    % Pk in PHedge
    PH_k = Pkappa(:, 2:2:end, :);
    PH_k(sign_H<0) = hardentanh2(PH_k(sign_H<0), sigma, PHedge1);
    PH_k(sign_H>0) = hardentanh2(PH_k(sign_H>0), sigma, PHedge2);
    PH_k(sign_H==0) = 0;
    Plk(:, 2:2:end, :, 1) = Plk(:, 2:2:end, :, 2) + PH_k.*(1-lambda);
    Plk(:, 2:2:end, :, 2) = PH_k;

    % reshape
    Plk = reshape(Plk, [], 2);
else
    Nslice2 = Nslice/2;
    Plk = cat(2, reshape(dataflow.rawdata(:, 1:2:end, :), [], 1), reshape(dataflow.rawdata(:, 2:2:end, :), [], 1));
    % 1/1000
    Plk = Plk./HU;
    % L/H to lamda/kappa
    Plk = Plk * Mlambda';
end
Np = size(Plk, 1);

% renorm Pl
Plk_sign = sign(Plk);
Plk(:, 1) = tablerenorm(Plk(:, 1), mdcorr.Plrange, mdcorr.NPlsamp);
% Plk(:, 1) = (Plk(:, 1) - mdcorr.Plrange(1)) ./ (mdcorr.Plrange(2) - mdcorr.Plrange(1));
% Plk(:, 1) = abs(Plk(:, 1)).*(mdcorr.NPlsamp - 1) + 1;
% Plk(Plk(:, 1) > mdcorr.NPlsamp, 1) = mdcorr.NPlsamp;

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
R = reshape(R, Npixel, Nslice2, Nview);
for islice = 1:Nslice2
    R_ii = squeeze(R(:, islice, :));
    R(:, islice, :) = BregmanTV2D(R_ii, 0.5, 0.1);
end

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
D = interp2(mdcorr.Tablelr2D, Rnorm(:), Plk(:,1));
D = D.*Plk_sign(:, 1);

% return rawdata
dataflow.rawdata = reshape(D.*(1+R(:).*1i).*HU, Npixel*Nslice2, Nview);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function Pout = tablerenorm(Pin, range, Nsamp)

Pout = (Pin - range(1)) ./ (range(2) - range(1));
Pout = abs(Pout).*(Nsamp - 1) + 1;
Pout(Pout > Nsamp) = Nsamp;

end
