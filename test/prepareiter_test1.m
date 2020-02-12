function [b0, b_noise] = prepare(x0, theta, emean)

b0 = radon(x0, theta);

mu = 0.05/1000;
In = 10000;
I = exp(-b0.*mu).*In;

b_noise = -log(poissrnd(I./emean).*emean./In)./mu;

