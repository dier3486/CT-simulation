function [dataflow, prmflow, status] = reconnode_materialdecompcali(dataflow, prmflow, status)
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

% parameters set in pipe
caliprm = prmflow.pipe.(status.nodename);

% % test manual settings
% caliprm.calidata = 'E:\matlab\CT\WY\CT3.0\calibration\tablesample\MDcalidata_newcorr.mat';
% caliprm.Plamdasamplenumber = 61;
% caliprm.Pkappasamplenumber = [3 12];
% caliprm.Rsamplenumber = [3 12];
% caliprm.toplot = true;

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

% electronic density or equavalent water density
if isfield(caliprm, 'elecdensity_flag')
    elecdensity_flag = caliprm.elecdensity_flag;
else
    elecdensity_flag = true;
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

% material physics
% effective Z and mass
[ZA, WAeff] = effectZandMass(materialA);
[ZB, WBeff] = effectZandMass(materialB);
% effective density
densityA = materialA.density;
densityB = materialB.density;
sigmaMZ = 2.0;      % norminal nuclear mass number / charge number
% sigmaMZ = WAeff/ZA;
if elecdensity_flag
    % effective density normalized mu which will be almost in Z^3 scaling.
    densityA = densityA * (ZA*sigmaMZ/WAeff);
    densityB = densityB * (ZB*sigmaMZ/WBeff);
end

% normalized mu
muA = materialA.mu_total;
muB = materialB.mu_total./densityB.*densityA;
% muAB
muAB = [muA  muB-muA];

% material C is used for calibration
if isfield(caliprm, 'calimaterial')
    materialC_cfg = fullfile(materialpath, caliprm.calimaterial);
else
    materialC_cfg = fullfile(materialpath, 'metalAl');
end
materialC = materialdefine(loadmaterial(materialC_cfg), samplekeV);
% effective ZC and density
[ZC, WCeff] = effectZandMass(materialC);
densityC = materialC.density;
if elecdensity_flag
    densityC = densityC * (ZC*sigmaMZ/WCeff);
end

% we do not use the mu_total of materialC
Rzc = (ZC^3 - ZA^3) / (ZB^3 - ZA^3);
% muZC = muA + (muB-muA).*Rzc;

% densitynormbyC = (densityC/densityCeff - densityA/densityAeff*(1-Rzc))/Rzc;
% I know the normalized dansity of B is densityBeff * densitynormbyC, which will be obviously large than densityB

% source sprectrum
% plz set one KV once
spectrum0 = prmflow.SYS.source.spectrum{1};

% lambda
if isfield(caliprm, 'LHmixlambda')
    LHmixlambda = caliprm.LHmixlambda;
else
    LHmixlambda = 0.3;
end
Mlambda = [1-LHmixlambda LHmixlambda; 1 -1];
% Minvlam = [1 LHmixlambda; 1 LHmixlambda-1];

% step0 
% the Beam harden correction based on material A
Dbh = -20:2:800;
% Pbh_L = -log((spectrum0.*resp_L.*samplekeV) * exp(-materialA.mu_total*Dbh) ./ sum(spectrum0.*resp_L.*samplekeV));
% Pbh_H = -log((spectrum0.*resp_H.*samplekeV) * exp(-materialA.mu_total*Dbh) ./ sum(spectrum0.*resp_H.*samplekeV));
Pbh = projectD(Dbh, samplekeV, spectrum0, resp, materialA.mu_total);

% V_lambda and V_kappa
Vab = (spectrum0.*samplekeV.*resp) * muAB;
Vlk = Mlambda*(Vab(:,2)./Vab(:,1));

% step1 calibration by experiment of material C
% calibration data
calidata = loaddata(caliprm.calidata);
% skip
calidata.D = calidata.D(1:end-1);
calidata.low = calidata.low(1:end-1);
calidata.high = calidata.high(1:end-1);

% experiment data
[~, snot0] = sort(calidata.D);
snot0 = snot0(calidata.D(snot0)~=0);
PLexp = calidata.low(snot0);
PHexp = calidata.high(snot0);
D0 = calidata.D(snot0);
Dexp = calidata.D(snot0).*(densityC/densityA);

% l: lambda
% k: kappa
% Pl Pk of experiments
Plkexp = Mlambda * [PLexp; PHexp];
% Plkexp(1, :) = calidata.low.*(1-LHmixlambda) + calidata.high.*LHmixlambda;
% Plkexp(2,:) = calidata.low - calidata.high;
Nexp = length(Dexp);

% fit R to Rzc
% by solving DcRfit = x: P(x, Rzc) = [Pl; *].
% and PkRfit = P(DcRfit, Rzc).
PkRfit = zeros(1, Nexp);
DcRfit = zeros(1, Nexp);
for ii = 1:Nexp
    DcRfit(ii) = fzero(@(x) projectPl(x, Rzc, samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) - Plkexp(1, ii), Plkexp(1, ii));
    PkRfit(ii) = projectPk(DcRfit(ii), Rzc, samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh);
end

% by mapping Plkexp(2,:) to PkRfit
% PbBound = max(Plkexp(1, :));
PlrBound = 120;
sigmaLR = 5e-3;
t1 = hardentanh2(Plkexp(1, :), sigmaLR, PlrBound);
g2r = polyfit(t1, PkRfit./Plkexp(2,:), 2);
% then the P(*, Rzc) ~= [Plkexp(1, :), Plkexp(2,:).*g2r(Plkexp(1, :))], that can be used to set up a table in looking for Rz by [Pl, Pk].

% fit D to Dexp
PldBound = 300;
sigmaLD = 5e-3;
Pldshift = 10;
t2 = log2(hardentanh2(Plkexp(1, :) + Pldshift, sigmaLD, PldBound));
g2d = polyfit(t2, DcRfit./Dexp-1, 3);

if toplot
    % plot step1 g2r
    tt = 1:500;
    figure;
    subplot(2,1,1); hold on
    plot(Plkexp(1, :), Plkexp(2,:), '.-');
    plot(Plkexp(1, :), PkRfit);
    grid on
    subplot(2,1,2); hold on
    plot(Plkexp(1, :), PkRfit./Plkexp(2,:), '.-');
    plot(tt, polyval(g2r, hardentanh2(tt, sigmaLR, PlrBound)));
    grid on;

    % plot step1 g2d
    tt = 1:500;
    figure;
    subplot(2,1,1); hold on
    plot(Plkexp(1, :), Dexp, '.-');
    plot(Plkexp(1, :), DcRfit);
    grid on
    subplot(2,1,2); hold on
    plot(Plkexp(1, :), DcRfit./Dexp, '.-');
    plot(tt, polyval(g2d, log2(hardentanh2(tt + Pldshift, sigmaLD, PldBound))) + 1);
    grid on
end

% step2 range of Pk, and table size
% Z range
if isfield(caliprm, 'Zrange')
    Zrange = caliprm.Zrange;
else
    % [1 50] is defualt range
    Zrange = [1 50];
end
% semi rate of Z 
Rzrange = (Zrange.^3 - ZA^3)./(ZB^3-ZA^3);
% muZrange = muA + (muB-muA).*Rzrange;

% Pl sample
if isfield(caliprm, 'Plamdarange')
    Plrange = caliprm.Plamdarange;
else
    Plrange = [0, 600];
end
if isfield(caliprm, 'Plamdasamplenumber')
    NPlsamp = caliprm.Plamdasamplenumber;
    Plsamp = linspace(Plrange(1), Plrange(2), NPlsamp)';
else
    Plsamp = (Plrange(1) : 2 : Plrange(2))';
    NPlsamp = length(Plsamp);
end

% Pkrange
Pkrange = zeros(NPlsamp, 2);
for ii = 1 : NPlsamp
    Dz1 = fzero(@(x) projectPl(x, Rzrange(1), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) - Plsamp(ii), Plsamp(ii));
    Dz2 = fzero(@(x) projectPl(x, Rzrange(2), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) - Plsamp(ii), Plsamp(ii)/Rzrange(2));
    Pkrange(ii, 1) = projectPk(Dz1, Rzrange(1), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh);
    Pkrange(ii, 2) = projectPk(Dz2, Rzrange(2), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh);
end
1;
Pkrange = Pkrange./polyval(g2r, hardentanh2(Plsamp, sigmaLR, PlrBound));

% Pk sample
if isfield(caliprm, 'Pkappasamplenumber')
    [NPksamp_neg, NPksamp_pos]= tac(caliprm.Pkappasamplenumber);
else
    NPksamp_neg = 21;
    NPksamp_pos = 300;
end
NPksamp = NPksamp_neg + NPksamp_pos - 1; 
Pksamp = zeros(NPlsamp, NPksamp);
Pksamp(:, 1:NPksamp_neg) = Pkrange(:, 1) * linspace(1, 0, NPksamp_neg);
Pksamp(:, NPksamp_neg : NPksamp) = Pkrange(:, 2) * linspace(0, 1, NPksamp_pos);

% step3 set up table [Pl; Pk] -> r
Tlk2r = zeros(NPlsamp, NPksamp);
Dlk = zeros(NPlsamp, NPksamp);
Pk2 = zeros(NPlsamp, NPksamp);
% solve r = x: P(*, x) = [Pl; Pk']. Pk' = Pk*g2r(Pl);
options = optimoptions('fsolve','Display','off');
% options = [];
% tic;
for ii =  1 : NPlsamp
    for jj = 1 : NPksamp
        Pl_ij = Plsamp(ii);
        Pk_ij = Pksamp(ii, jj)*polyval(g2r, hardentanh2(Pl_ij, sigmaLR, PlrBound));
        if ii == 1
            u0 = [0 0];
        else
            u0 = [Dlk(ii-1, jj) Tlk2r(ii-1, jj)];
        end
        [u, ~, fsflag] = fsolve(@(x) projectP(x(1), x(2), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) ...
            - [Pl_ij; Pk_ij], u0, options);
        if ~fsflag
            1;
        end
        Dlk(ii, jj) = u(1);
        Tlk2r(ii, jj) = u(2);
        
        Pk2(ii, jj) = Pk_ij;
    end
end
% toc;
% 
% if s0 > 1
%     Tlk2r(s0, :) = (Tlk2r(s0-1, :) + Tlk2r(s0+1, :))./2;
% end

% when Plsamp=0,
s0 = find(Plsamp == 0);
Tlk2r(s0, 1 : NPksamp_neg) = 1./(1./(NPksamp_neg-1:-1:0).*(NPksamp_neg-1).*( Vlk(1) + 1/Rzrange(1)) - Vlk(1));
Tlk2r(s0, NPksamp_neg : NPksamp) = 1./(1./(0:NPksamp_pos-1).*(NPksamp_pos-1).*( Vlk(1) + 1/Rzrange(2)) - Vlk(1));

% plot step3 Tlk2r
if toplot
    figure;
    mesh(repmat(Plsamp, 1, NPksamp), Pksamp, Tlk2r);
    grid on;
end

% step4 setup table [Pl; r] -> Drho
% table size
if isfield(caliprm, 'Rsamplenumber')
    [NRsamp_neg, NRsamp_pos]= tac(caliprm.Rsamplenumber);
else
    NRsamp_neg = 21;
    NRsamp_pos = 300;
end
NRsamp = NRsamp_neg + NRsamp_pos - 1;
% R sampling
Rsamp = zeros(1, NRsamp);
Rsamp(1:NRsamp_neg) = Rzrange(1) .* linspace(1, 0, NRsamp_neg);
Rsamp(NRsamp_neg : NRsamp) = Rzrange(2) .* linspace(0, 1, NRsamp_pos);
% Zsamp = (ZA^3.*(1-Rsamp) + ZB^3.*Rsamp).^(1/3);
Tlr2D = zeros(NPlsamp, NRsamp);
Dlr = zeros(NPlsamp, NRsamp);
% solve D = x: P(x, r) = [Pl; *].
Dscale_R = (1/Rzc + Vlk(1))./(1./Rsamp + Vlk(1));
% tic;
for ii = 1 : NPlsamp
    for jj = 1:NRsamp
        if ii == 1
            u0 = 0;
        else
            u0 = Dlr(ii-1, jj);
        end
        Dlr(ii, jj) = fzero(@(x) projectPl(x, Rsamp(jj), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) ...
            -Plsamp(ii), u0);
        
%         D0 = Dlr(ii, jj) / (1 + polyval(g2d, log2(hardentanh2(Plsamp(ii) + Pldshift, sigmaLD, PldBound))) * (Zsamp(jj)-ZA) / (ZC-ZA) );
        D0 = Dlr(ii, jj) / (1 + polyval(g2d, log2(hardentanh2(Plsamp(ii) + Pldshift, sigmaLD, PldBound))) * Dscale_R(jj) );
        Tlr2D(ii, jj) = D0;
    end
end
% toc;

% plot step4 Tlr2D
if toplot
    figure;
    mesh(repmat(Plsamp, 1, NRsamp), repmat(Rsamp, NPlsamp, 1), Tlr2D);
    grid on;
end

% step5 edge
if isfield(caliprm, 'Zedge')
    Zedge = caliprm.Zedge;
else
    Zedge = Zrange;
end
% semi rate of Z 
Rzedge = (Zedge.^3 - ZA^3)./(ZB^3-ZA^3);
% Pkedge
Pkedge = zeros(NPlsamp, 2);
for ii = 1 : NPlsamp
    Dz1 = fzero(@(x) projectPl(x, Rzedge(1), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) - Plsamp(ii), Plsamp(ii));
    Dz2 = fzero(@(x) projectPl(x, Rzedge(2), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh) - Plsamp(ii), Plsamp(ii)/Rzedge(2));
    Pkedge(ii, 1) = projectPk(Dz1, Rzedge(1), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh);
    Pkedge(ii, 2) = projectPk(Dz2, Rzedge(2), samplekeV, spectrum0, resp, muAB, LHmixlambda, Dbh, Pbh);
end
Pkedge = Pkedge./polyval(g2r, hardentanh2(Plsamp, sigmaLR, PlrBound));
% PL edge, when PL is fixed, the range of Pk
PLedge = zeros(NPlsamp, 2);
for ii = 1 : NPlsamp
    if ii == 1
        x0 = 0;
    else
        x0 = PLedge(ii-1, 1);
    end
    PLedge(ii, 1) = fzero(@(x) interp1(Plsamp, Pkedge(:, 1), Plsamp(ii) - x*LHmixlambda, ...
        'linear', 'extrap') - x, x0);
    if ii == 1
        x0 = 0;
    else
        x0 = PLedge(ii-1, 2);
    end
    PLedge(ii, 2) = fzero(@(x) interp1(Plsamp, Pkedge(:, 2), Plsamp(ii) - x*LHmixlambda, ...
        'linear', 'extrap') - x, x0);
end


% PH edge, , when PH is fixed, the range of Pk
PHedge = zeros(NPlsamp, 2);
for ii = 1 : NPlsamp
    if ii == 1
        x0 = 0;
    else
        x0 = PHedge(ii-1, 1);
    end
    PHedge(ii, 1) = fzero(@(x) interp1(Plsamp, Pkedge(:, 1), Plsamp(ii) + x*(1-LHmixlambda), ...
        'linear', 'extrap') - x, x0);
    if ii == 1
        x0 = 0;
    else
        x0 = PHedge(ii-1, 2);
    end
    PHedge(ii, 2) = fzero(@(x) interp1(Plsamp, Pkedge(:, 2), Plsamp(ii) + x*(1-LHmixlambda), ...
        'linear', 'extrap') - x, x0);
end

% calibration table
mdcorr = caliprmforcorr(prmflow, corrversion);
mdcorr.elecdensity_flag = elecdensity_flag;
mdcorr.ZA = ZA;
mdcorr.densityA = densityA;
mdcorr.ZB = ZB;
mdcorr.densityB = densityB;
mdcorr.ZC = ZC;
mdcorr.densityC = densityC;
mdcorr.LHmixlambda = LHmixlambda;
mdcorr.Zrange = Zrange;
mdcorr.Rrange = Rzrange;
mdcorr.Plrange = Plrange;
mdcorr.NPlsamp = NPlsamp;
mdcorr.Pkrange = Pkrange;  % in size NPlsamp*2
mdcorr.NPksamp = [NPksamp_neg NPksamp_pos];
mdcorr.NRsamp = [NRsamp_neg  NRsamp_pos];
mdcorr.Tablelk2r = Tlk2r;  % in size NPlsamp * (NPksamp(1) + NPksamp(2) - 1)
mdcorr.Tablelr2D = Tlr2D;  % in size NPlsamp * (NRsamp(1) + NRsamp(2) - 1)
mdcorr.Zedge= Zedge;
mdcorr.PLedge = PLedge;  % in size NPlsamp*2
mdcorr.PHedge = PHedge;  % in size NPlsamp*2
mdcorr.Vlk = Vlk;


% to single
mdcorr = everything2single(mdcorr, 'double', 'single');

% to return
dataflow.materialdecompcorr = mdcorr;

end

% sub functions
function P = projectD(Drho, samplekeV, spectrum0, resp, mu)

% I know resp = [resp_L; resp_H];

P = -log((spectrum0.*resp.*samplekeV) * exp(-mu*Drho)) + log(sum(spectrum0.*resp.*samplekeV, 2));
P(P<-100) = -100; 

end

function Pl = projectPl(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
Pl = [1-lambda  lambda] * P2;

end

function Pk = projectPk(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
Pk = [1  -1] * P2;

end

function P = projectP(Drho, r, samplekeV, spectrum0, resp, mu, lambda, Dbh, Pbh)

mu = mu * [1; r];
P1 = projectD(Drho, samplekeV, spectrum0, resp, mu);
P2 = [interp1(Pbh(1, :), Dbh, P1(1, :), 'linear', 'extrap'); 
      interp1(Pbh(2, :), Dbh, P1(2, :), 'linear', 'extrap')];
P = [1-lambda  lambda; 1  -1] * P2;

end