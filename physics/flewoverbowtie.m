function [Dmu, L] = flewoverbowtie(A, B, bowtie, filter, samplekeV)
% projection on bowtie and filter, from position A to position B
% Dmu = flewoverbowtie(A, B, bowtie, filter, samplekeV);
% or typically,
% [Dmu, L] = flewoverbowtie(focalposition, detectorposition, bowtie, filter, samplekeV);
% NOTE: for B in (m,3) and A in (n,3) the returnd L is in (m,n) and Dmu is in (m*n, Nsample)
% where the n is, typically, the focal number.

% x y z
xx = B(:,1) - A(:,1)';
yy = B(:,2) - A(:,2)';
zz = B(:,3) - A(:,3)';
nA = size(A, 1);
% tube number and focal number
% geometry
XYangle = atan2(yy, xx) - atan2(-A(:,2), -A(:,1))';
XYangle = mod(XYangle+pi, pi*2)-pi;
% Zscale = sqrt(yy.^2+zz.^2)./yy; 
Zscale = sqrt(xx.^2 + yy.^2+zz.^2)./sqrt(xx.^2 + yy.^2);
% Dfscale = (sqrt(xx.^2+yy.^2+zz.^2)./yy);

% ini Dmu
Nd = size(xx(:), 1);
Nsample = length(samplekeV(:));
Dmu = zeros(Nd, Nsample);

% bowtie(s)
Nbowtie = length(bowtie);
for ibow = 1:Nbowtie
    bowtie_ii = bowtie{ibow};
    % I know the focal spots number is
    Nfspot = size(bowtie_ii.bowtiecurve, 2);
    if isempty(bowtie_ii.bowtiecurve)
        % empty bowtie
        continue;
    end
    % D
    D_bowtie = zeros(size(XYangle));
    for jj = 1:nA
        ispot = mod(jj-1, Nfspot)+1;
        D_bowtie(:, jj) = interp1(bowtie_ii.anglesample(:, ispot), double(bowtie_ii.bowtiecurve(:, ispot)), ...
            XYangle(:, jj), 'linear', 'extrap');
    end
    D_bowtie = D_bowtie.*Zscale;
    % mu
    mu_bowtie = interp1(bowtie_ii.material.samplekeV, bowtie_ii.material.mu_total, samplekeV);
    % I know in most case bowtie.material.samplekeV == samplekeV
    % + to Dmu
    Dmu = Dmu + D_bowtie(:)*mu_bowtie;
end

% filter(s)
Nfilter = length(filter(:));
for ifil = 1:Nfilter
    filter_ii = filter{ifil};
    if ~isfield(filter_ii, 'origangle')
        filter_ii.origangle = zeros(Nfspot,1);
    end
    % D
    if isfield(filter_ii, 'effect') && filter_ii.effect
        % do not scale by angle;
        D_filter = filter_ii.thickness;
        % I know this line will be used in BH cali.
    else
        for jj = 1:nA
            ispot = mod(jj-1, Nfspot)+1;
        end
        Dfscale = sec(XYangle(:, jj) + filter_ii.origangle(ispot)).*Zscale;
        D_filter = Dfscale.*filter_ii.thickness;
    end
    % mu
    mu_filter = interp1(filter_ii.material.samplekeV, filter_ii.material.mu_total, samplekeV);
    % + to Dmu
    Dmu = Dmu + D_filter(:)*mu_filter;
end

% drop by L 
L = sqrt(xx.^2+yy.^2+zz.^2);

end