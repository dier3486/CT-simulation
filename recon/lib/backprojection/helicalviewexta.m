function [viewextra_full, viewextra_pi] = helicalviewexta(Nslice, pitch, effRadius, Nviewprot)
% subfunction of helicalprepare.m and cradleprepare.m

sigma_z = (Nslice-1)/Nslice;
absPitch = abs(pitch);

t0 = fzero(@(x) (effRadius + sqrt((1-effRadius^2*x^2)/(1-x^2)) )*x*sigma_z - absPitch/pi, 0);
Z0 = asin(effRadius*t0) * absPitch/(pi*2) + (sqrt(1-effRadius^2*t0^2) + effRadius*sqrt(1-t0^2)) * sigma_z/2;
viewextra_full = Z0/absPitch*Nviewprot;

% pi-line recon condition
t_pi = fzero(@(x)  (effRadius^2 * x * sqrt(1-x^2)/(1-effRadius^2*x^2) - x/sqrt(1-x^2))*pi/2 - 1, 0);
Zpi = (1+effRadius*sqrt( (1-t_pi^2)/(1-effRadius^2*t_pi^2) ) )/4 - asin(t_pi*effRadius)/pi/2;
viewextra_pi = Zpi*Nviewprot;

end