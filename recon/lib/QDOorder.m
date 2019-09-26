function [a1, a2] = QDOorder(N, m)

% % debug
% N = 11;
% m = 4.13;
% dd = (1:N).' - m;

minl = floor((1-m)*2);
minr = floor((m-N)*2);
maxl = floor((N-m)*2);
maxr = floor((m-1)*2);

minlr = max(minl, minr) - 1;
maxlr = min(maxl, maxr) + 1;

xx = 1:N;
a1 = floor((xx-m).*2);
a2 = floor((-xx+m).*2);
a1(a1<minlr | a1>maxlr) = nan;
a2(a2<minlr | a2>maxlr) = nan;
a1 = a1 - minlr + 1;
a2 = a2 - minlr + 1;

end

% % debug
% s1 = ~isnan(a1);
% s2 = ~isnan(a2);
% N_QDO = max([a1, a2]);
% d_QDO = nan(N_QDO, 1);
% d_QDO(a1(s1)) = dd(s1);
% d_QDO(a2(s2)) = -dd(s2);

