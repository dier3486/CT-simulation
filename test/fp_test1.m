Cimage = phantom().*100;

[Nx, Ny] = size(Cimage);

Np = 500;
A = [0 -200 0];
B = [linspace(-200, 200, Np) ones(Np,1).*200 zeros(Np,1)];
h = 1;

theta = 0.1;