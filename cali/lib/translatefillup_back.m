function [raw1, n_l] = translatefillup(raw0, len1, mid_u, blkvindex)
% fillup rawdata by translate the center data,
% raw1 = translatefillup(raw0, len1, mid_u, blkvindex);
% or, forced to fill up all the views
% raw1 = translatefillup(raw0, len1, mid_u, true(1, Nview));

[Npixel, Nview] = size(raw0);
if nargin<3
    mid_u = Npixel/2;
end
mid_u = double(mid_u);
n_l = floor((len1-Npixel)/2);
n_r = ceil((len1-Npixel)/2);

% ini
raw1 = zeros(len1, Nview);
raw1(n_l+1:len1-n_r, :) = raw0;

if ~any(blkvindex)
    return
end

% find a center projection
x0 = 1:Npixel;
w1 = double(x0*raw0./sum(raw0, 1));
[~, s_mid] = min(abs(w1-mid_u));
a0 = raw0(:,s_mid);

% fillup
x_fl = [(1:n_l)-n_l  (1:n_r)+Npixel];
s_fl = x_fl+n_l;
viewindex = 1:Nview;
options_fzero = optimset('TolX',1e-8);
for ii = viewindex(blkvindex)
    x = fzero(@(t) alignfit(a0, w1(ii), t), w1(ii)-mid_u, options_fzero);
    raw1(s_fl, ii) = interp1(x0, a0, x_fl-x,'linear', 0);
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