function y = mylorenzcurve(x, alpha)
% a Lorenz curve, hyperbolic
% y = mylorenzcurve(x, alpha);
% x in [0 1], alpha in [-1, 1].

y = (1-alpha).*x./(1+alpha - x.*alpha.*2);