function M = backproj2D_spM1(Np, theta, ctrIdx, d, N)
% back projection sparse matrix (for linear interp method)
% M = backproj2D_1(Np, theta, ctrIdx, d, N, interp)

% Define the x & y axes for the reconstructed image so that the origin
% (center) is in the spot which RADON would choose.
center = floor((N + 1)/2);
xleft = -center + 1;
x = (1:N) - 1 + xleft;
x = repmat(x, N, 1);

ytop = center - 1;
y = (N:-1:1).' - N + ytop;
y = repmat(y, 1, N);

% Generate trignometric tables
costheta = cos(theta);
sintheta = sin(theta);

Nview = length(theta);
M = sparse([]);
% interp
index_p = 1:N^2;
for iview=1:length(theta)
    t = (x.*costheta(iview) + y.*sintheta(iview)).*d + ctrIdx;
    t_index1 = floor(t);
    t_alpha = t - t_index1;
    t_index2 = t_index1+1;
    s1 = t_index1>0 & t_index1<=Np;
    s2 = t_index2>0 & t_index2<=Np;
    M1 = sparse(index_p(s1), t_index1(s1), 1-t_alpha(s1), N^2, Np);
    M2 = sparse(index_p(s2), t_index2(s2), t_alpha(s2), N^2, Np);
    M = [M, M1+M2];
end
        
M = M.*(pi/(2*Nview));

return