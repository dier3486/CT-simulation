% load E:\matlab\CT\SINO\PG\smiletest1.mat

tshell = 5;
tscale = 1.2;
delta_d = double(tshell./prmflow.recon.delta_d*(tscale-1));

n1 = find(Amean>C0, 1, 'first');
n2 = find(Amean>C0, 1, 'last');

d1 = n1 - (Amean(n1) - C0)/(Amean(n1) - Amean(n1-1)) - delta_d;
d2 = n2 + (Amean(n2) - C0)/(Amean(n2) - Amean(n2+1)) + delta_d;

d0 = (d1+d2)/2;
r0 = (d2-d1)/2;

m = 0;
xx = find(Amean(n1+m : n2-m) < C0) + n1 + m - 1;
tt = [d1 n1:n2 d2];

% a = sum(C0-Amean(xx))/sum(sqrt(r0^2 - (xx-d0).^2));
% yt = C0 - sqrt(r0^2 - (tt-d0).^2).*a;

% a = lsqnonlin(@(a) (cosh((xx-d0)./r0./a) - cosh(1/a)).*a - (Amean(xx)-C0), 0.1);
% yt = C0 + (cosh((tt-d0)./r0./a) - cosh(1/a)).*a;

% a = sum(Amean(xx)-C0)/sum( (1+(xx-d0)./r0).*(1-(xx-d0)./r0));
% yt = C0 + (1+(tt-d0)./r0).*(1-(tt-d0)./r0).*a;
options = optimoptions('lsqnonlin','Display','off');
a = lsqnonlin(@(a) (1+(xx-d0)./r0).*(1-(xx-d0)./r0).*(a(1)+a(2).*((xx-d0)./r0).^2) - (Amean(xx)-C0), [0.1 0], [], [], options);
yt = C0 + (1+(tt-d0)./r0).*(1-(tt-d0)./r0).*(a(1)+a(2).*((tt-d0)./r0).^2);

figure;
hold on
plot(Amean);
plot([d1 d2], [C0 C0], '*-');
plot(tt, yt);
% axis([n1-10, n2+10, C0-50 C0+20]);
grid on