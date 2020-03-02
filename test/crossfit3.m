function r = crossfit3(p, x, y, s, lambda)

r1 = x(2,:) + p(1).*y(1,:) + p(2).*y(3,:) - y(2,:).*(1+p(1)+p(2));
r1 = r1./(-y(2,:).*log2(y(2,:)));

r = [r1(s)  p.*lambda];

end