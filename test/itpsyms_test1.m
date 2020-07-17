% syms x_2 x_1 x0 x1 x2
% syms a c
% 
% b = 1/2-1/2*sqrt(1+4*a*(1-a));
% y0 = (1-a+b)*x_1/4 + (2-a-b)*x0/4 + (1+a-b)*x1/4 + (a+b)*x2/4;
% au = a+1/2;
% bu = 1/2-1/2*sqrt(1+4*au*(1-au));
% yu = (1-au+bu)*x_1/4 + (2-au-bu)*x0/4 + (1+au-bu)*x1/4 + (au+bu)*x2/4;
% yd = (1-au+bu)*x_2/4 + (2-au-bu)*x_1/4 + (1+au-bu)*x0/4 + (au+bu)*x1/4;
% 
% z = -yd*c/2 + y0*(1-c) -yu*c/2;

g1 = 0:0.01:1.5;
g2 = zeros(size(g1));
for ii = 1:length(g1)
    g2(ii) = fzero(@(x) (g1(ii)/8)^2+(1/4-x/4)^2+(1/2+g1(ii)/4+x*3/8)^2+1/16+(g1(ii)/8+x/8)^2-1, 0.5);
end