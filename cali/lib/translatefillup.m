function [raw1, n_l] = translatefillup(raw0, len1, mid_u)
% fillup rawdata by translate the center data,
% raw1 = translatefillup(raw0, len1, mid_u);

[Npixel, Nview] = size(raw0);
if nargin<3
    mid_u = Npixel/2;
end
mid_u = double(mid_u);
n_l = floor((len1-Npixel)/2);
n_r = ceil((len1-Npixel)/2);

% find a center projection
x0 = 1:Npixel;
w1 = double(x0*raw0./sum(raw0, 1));
[~, s_mid] = min(abs(w1-mid_u));
a0 = raw0(:,s_mid);

% fillup
raw1 = zeros(len1, Nview);
raw1(n_l+1:len1-n_r, :) = raw0;
x_fl = [(1:n_l)-n_l  (1:n_r)+Npixel];
s_fl = x_fl+n_l;
for ii = 1:Nview
    x = fzero(@(t) alignfit(a0, w1(ii), t), w1(ii)-mid_u);
    raw1(s_fl, ii) = interp1(x0, a0, x_fl-x);
end
raw1(isnan(raw1)) = 0;

end