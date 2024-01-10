% function [dataflow, prmflow, status] = reconnode_materialdecompcali(dataflow, prmflow, status)
% recon node, beamharden calibration
% [dataflow, prmflow, status] = reconnode_materialdecompcali(dataflow, prmflow, status);

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

load F:\物影\MaterialDecomp\decomp1.mat

% parameters set in pipe
caliprm = prmflow.pipe.(status.nodename);

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end
% to plot 
if isfield(caliprm, 'toplot')
    toplot = caliprm.toplot;
else
    toplot = false;
end
% resp
samplekeV = prmflow.SYS.world.samplekeV;
if isfield(caliprm, 'responseLow')
    respdata = load(caliprm.responseLow);
    resp_L = interp1(respdata.samplekeV, respdata.response, samplekeV);
else
    % I know
    resp_L = prmflow.SYS.detector.response(1, :);
end
if isfield(caliprm, 'responseHigh')
    respdata = load(caliprm.responseHigh);
    resp_H = interp1(respdata.samplekeV, respdata.response, samplekeV);
else
    % I know
    resp_H = prmflow.SYS.detector.response(end, :);
end
resp = [resp_L; resp_H];

% % debug
% metalCu = materialdefine(loadmaterial('metalCu'), samplekeV);
% resp = [ones(size(samplekeV)); exp(-metalCu.mu_total.*1.5)'];
resp(resp<eps) = eps;

% base material A B
materialpath = prmflow.SYS.world.materialdata;
if isfield(caliprm, 'basematerialA')
    materialA_cfg = fullfile(materialpath, caliprm.basematerialA);
else
    % water is base material A
    materialA_cfg = fullfile(materialpath, 'water');
end
materialA = materialdefine(loadmaterial(materialA_cfg), samplekeV);

if isfield(caliprm, 'basematerialB')
    materialB_cfg = fullfile(materialpath, caliprm.basematerialB);
else
    % the default base material B is Titanium
    materialB_cfg = fullfile(materialpath, 'metalTi');
end
materialB = materialdefine(loadmaterial(materialB_cfg), samplekeV);
% Z of the base materials
if isfield(caliprm, 'baseZA')
    ZA = caliprm.baseZA;
    % plz set water to 7.42;
else
    ZA = (([materialA.elemprm(:).Z].^4*materialA.elemmol) / ([materialA.elemprm(:).Z]*materialA.elemmol)).^(1/3);
    % ZA(water) = 7.4278;
end
if isfield(caliprm, 'baseZB')
    ZB = caliprm.baseZB;
else
    ZB = (([materialB.elemprm(:).Z].^4*materialB.elemmol) / ([materialB.elemprm(:).Z]*materialB.elemmol)).^(1/3);
end

% material physics (for test)
% sumZA = 
sumWA = [materialA.elemprm(:).weight] * materialA.elemmol;
sumWB = [materialB.elemprm(:).weight] * materialB.elemmol;
[ZAeff, WAeff] = effectZandMass(materialA);
[ZBeff, WBeff] = effectZandMass(materialB);
densityAeff = materialA.density * (ZAeff*2/WAeff);
densityBeff = materialB.density * (ZBeff*2/WBeff);

% Z sample
if isfield(caliprm, 'Zrange')
    Zrange = caliprm.Zrange;
else
    % [1 50] is defualt range
    Zrange = [1 50];
end
Zstep = 5;
Zsamp = Zrange(1) : Zstep : Zrange(2);
Nz = length(Zsamp);

% plz set one KV once
spectrum0 = prmflow.SYS.source.spectrum{1};

% density normalized mu
muA = materialA.mu_total;
muB = materialB.mu_total.*(materialA.density/materialB.density);
muAB = [muA  muB-muA];

% semi rate of Z 
Rzrange = (Zrange.^3 - ZA^3)./(ZB^3-ZA^3);
muZrange = muA + (muB-muA).*Rzrange;

% Beam harden correction based on material A
Dbh = 0:2:600;
% Pbh_L = -log((spectrum0.*resp_L.*samplekeV) * exp(-materialA.mu_total*Dbh) ./ sum(spectrum0.*resp_L.*samplekeV));
% Pbh_H = -log((spectrum0.*resp_H.*samplekeV) * exp(-materialA.mu_total*Dbh) ./ sum(spectrum0.*resp_H.*samplekeV));
Pbh = projectD(Dbh, samplekeV, spectrum0, resp, materialA.mu_total);

% lambda
if isfield(caliprm, 'mixlambda')
    lambda = caliprm.mixlambda;
else
    lambda = 0.3;
end

NPbsamp = 10;
Pbsamp = linspace(0, 500, NPbsamp)';
Pnrange = zeros(NPbsamp, 2);
DBrange = zeros(NPbsamp, 2);

% step1 range of Pn
for ii = 2 : NPbsamp
    Dz1 = fzero(@(x) projectPb(x, Rzrange(1), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) - Pbsamp(ii), Pbsamp(ii));
    Dz2 = fzero(@(x) projectPb(x, Rzrange(2), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) - Pbsamp(ii), Pbsamp(ii)/Rzrange(2));
    Pnrange(ii, 1) = projectPn(Dz1, Rzrange(1), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh);
    Pnrange(ii, 2) = projectPn(Dz2, Rzrange(2), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh);
    DBrange(ii, 1) = Dz1 * Rzrange(1);
    DBrange(ii, 2) = Dz2 * Rzrange(2);
end

NPnsamp_neg = 5;
NPnsamp_pos = 10;
NPnsamp = NPnsamp_neg + NPnsamp_pos - 1;
Pnsamp = zeros(NPbsamp, NPnsamp);
Pnsamp(:, 1:NPnsamp_neg) = Pnrange(:, 1) * linspace(1, 0, NPnsamp_neg);
Pnsamp(:, NPnsamp_neg : NPnsamp) = Pnrange(:, 2) * linspace(0, 1, NPnsamp_pos);

% a21 = projectP(Pbsamp(2), 1-Pnsamp(2, 1)/Pbsamp(2), samplekeV, spectrum0, resp, [muA muB], lambda, Dbh, Pbh);
options = optimoptions('fsolve','Display','off');

% step2 setup table of f1g and f1c (or f1r)
f1b = zeros(NPbsamp, NPnsamp);    % f1c(Pb, Pn) = Db
f1r = zeros(NPbsamp, NPnsamp);    % f1r(Pb, Pn) = Rz
f1d = zeros(NPbsamp, NPnsamp);    % f1g(Pb, Pn) = D = Da + Db*rho_b/rho_a
for ii = 2 : NPbsamp
    for jj = 1 : NPnsamp
        u = fsolve(@(x) projectP(x(1), x(2), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) ...
            - [Pbsamp(ii); Pnsamp(ii, jj)], [Pbsamp(ii) 1-Pnsamp(ii, jj)/Pbsamp(ii)], options);
        f1d(ii, jj) = u(1);
        f1b(ii, jj) = u(1)*u(2);
        f1r(ii, jj) = u(2);
    end
end
% I know the P_B is f1b or f1d*f1r to be denoised
% after denoising the new r' is P_B'/f1d

% step4 setup table of f2d (or f2g)
f2d = zeros(NPbsamp, NPnsamp);    % f2d(Pb, Dc) = DT
f2g = zeros(NPbsamp, NPnsamp);    % f2g(Pb, Rz) = D
NRsamp_neg = 6;
NRsamp_pos = 9;
NRsamp = NRsamp_neg + NRsamp_pos - 1;
Rsamp = zeros(1, NRsamp);
Rsamp(1:NRsamp_neg) = Rzrange(1) .* linspace(1, 0, NRsamp_neg);
Rsamp(NRsamp_neg : NRsamp) = Rzrange(2) .* linspace(0, 1, NRsamp_pos);
DBsamp = zeros(NPbsamp, NRsamp);
DBsamp(:, 1:NRsamp_neg) = DBrange(:, 1) * linspace(1, 0, NRsamp_neg);
DBsamp(:, NRsamp_neg : NRsamp) = DBrange(:, 2) * linspace(0, 1, NRsamp_pos);

for ii = 2 : NPbsamp
    for jj = 1:NRsamp
        f2g(ii, jj) = fzero(@(x) projectPb(x, Rsamp(jj), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) ...
            -Pbsamp(ii), Pbsamp(ii)/(1+Rsamp(jj)));
        f2d(ii, jj) = fzero(@(x) projectPb(x, DBsamp(ii, jj)/x, samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) ...
            -Pbsamp(ii), Pbsamp(ii)/(1+Rsamp(jj)));
    end
end

% step4 cali by C
calidata = load('E:\matlab\CT\WY\CT3.0\calibration\tablesample\MDcalidata2.mat');
calidata.low = calidata.low(calidata.D~=0);
calidata.high = calidata.high(calidata.D~=0);
calidata.D = calidata.D(calidata.D~=0);

if isfield(caliprm, 'caliZC')
    ZC = caliprm.caliZC;
elseif isfield(caliprm, 'calimaterial')
    materialC_cfg = fullfile(materialpath, caliprm.calimaterial);
    materialC = materialdefine(loadmaterial(materialC_cfg), samplekeV);
    % ZC and density
    ZC = [materialC.elemprm(:).Z] * materialC.elemweight;
    densityC = materialC.density;
    % we do not use the mu_total of materialC
else
    % defualt cali material is Aluminium
    ZC = 13;
    densityC = 2.6941;
end
Rzc = (ZC^3 - ZA^3) / (ZB^3 - ZA^3);
muZC = muA + (muB-muA).*Rzc;
% NOTE: the muZC is not the mu_total of the cali material
PbCexp = calidata.low.*(1-lambda) + calidata.high.*lambda;
PnCexp = calidata.low - calidata.high;
Nexp = length(PbCexp);

% fit R to Rzc
PnRfit = zeros(1, Nexp);
for ii = 1:Nexp
    Dc = fzero(@(x) projectPb(x, Rzc, samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) - PbCexp(ii), PbCexp(ii));
    PnRfit(ii) = projectPn(Dc, Rzc, samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh);
end
% by mapping PnCexp to Pn
PbBound = max(PbCexp);
sigmaB = 1e-5;
t = hardentanh(PbCexp, sigmaB, PbBound);
g2r = polyfit(t, PnRfit./PnCexp, 2);

% fit D to Dexp
Dexp1 = calidata.D.*(densityC/materialA.density);
Dexp2 = calidata.D.*Rzc;

PnDfit = zeros(1, Nexp);
RDfit = zeros(1, Nexp);
for ii = 1:Nexp
    RDfit(ii) = fzero(@(x) projectPb(Dexp2(ii), x, samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh) - PbCexp(ii), Rzc);
    PnDfit(ii) = projectPn(Dexp(ii), RDfit(ii), samplekeV, spectrum0, resp, muAB, lambda, Dbh, Pbh);
end

% fit calidata.D.*densityC to PCfit(1, :).

Pprin1 = projectD(calidata.D, samplekeV, spectrum0, resp, materialC.mu_total);
Pprin2 = [interp1(Pbh(1, :), Dbh, Pprin1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, Pprin1(2, :), 'linear', 'extrap')];


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
% end


function P = projectD(Drho, samplekeV, spectrum0, resp, mu)

% I know resp = [resp_L; resp_H];

P = -log((spectrum0.*resp.*samplekeV) * exp(-mu*Drho)) + log(sum(spectrum0.*resp.*samplekeV, 2));
P(P<-100) = -100; 

end

function Pb = projectPb(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
Pb = [1-lambda  lambda] * P2;

end

function Pn = projectPn(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
Pn = [1  -1] * P2;

end

function P = projectP(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
P = [1-lambda  lambda; 1  -1] * P2;

end

