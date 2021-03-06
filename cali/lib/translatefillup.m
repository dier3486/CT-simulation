function [raw1, n_l] = translatefillup(raw0, len1, refblock)
% fillup rawdata by translate the center data,
% raw1 = translatefillup(raw0, len1, refblock);
% or, forced to fill up all the views
% raw1 = translatefillup(raw0, len1, mid_u, true(1, Nview));

[Npixel, Nview] = size(raw0);

n_l = floor((len1-Npixel)/2);
n_r = ceil((len1-Npixel)/2);

% ini
raw1 = zeros(len1, Nview);
raw1(n_l+1:len1-n_r, :) = raw0;

% nothing to fill?
if ~any(refblock(:))
    return;
end

% find a center projection
x0 = 1:Npixel;
wcenter = double(x0*raw0./sum(raw0, 1));

% index to fillup
x_left = (1:n_l)-n_l;
s_left = x_left + n_l;
x_right = (1:n_r)+Npixel;
s_right = x_right + n_l;

raw1 = fitandfill(raw1, raw0, refblock(1, :), wcenter, x_left, s_left);
raw1 = fitandfill(raw1, raw0, refblock(2, :), wcenter, x_right, s_right);

end

function A1 = fitandfill(A1, A0, refblock, wcenter, x_fl, s_fl)

[Npixel, Nview] = size(A0);
x0 = 1:Npixel;
options_fzero = optimset('TolX',1e-8);
for ii = find(refblock)
    v1 = find(~refblock(1:ii), 1, 'last');
    v2 = find(~refblock(ii+1:end), 1, 'first') + ii;
    if isempty(v1)
        v1 = find(~refblock, 1, 'last');
        d1 = ii + Nview - v1;
    else
        d1 = ii - v1;
    end
    if isempty(v2)
        v2 = find(~refblock, 1, 'first');
        d2 = v2 + Nview - ii;
    else
        d2 = v2 - ii;
    end
    % fit
    a1 = A0(:, v1);
    a2 = A0(:, v2);
    x1 = fzero(@(t) alignfit(a1, wcenter(ii), t), wcenter(ii)-wcenter(v1), options_fzero);
    x2 = fzero(@(t) alignfit(a2, wcenter(ii), t), wcenter(ii)-wcenter(v2), options_fzero);
    % fill
    A1(s_fl, ii) = interp1(x0, a1, x_fl-x1,'linear', 0).*(d2/(d1+d2)) + interp1(x0, a2, x_fl-x2,'linear', 0).*(d1/(d1+d2));
end

end

function r = alignfit(y, w, x)
% x = fzero(@(x) alignfit(y, w, x), x0);

y = y(:)';
N = length(y);
xx = 1:N;

y2 = interp1(xx, y, xx-x, 'linear', 0);
w2 = sum(y2.*xx)/sum(y2);
r = w2-w;

end