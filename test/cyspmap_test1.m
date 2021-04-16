tt = linspace(0, pi, 200);
f1 = zeros(1,200);
s1 = tt<=pi/4;
f1(s1) = tan(tt(s1)).^2./6;
s2 = tt>pi/4 & tt<=pi*3/4;
f1(s2) = (1-cot(tt(s2)))./3+1/6;
s3 = tt>pi*3/4;
f1(s3) = 1-tan(tt(s3)).^2./6;

xx = linspace(0,1,200);
f2 = zeros(1, 200);
s1 = xx<=1/6;
f2(s1) = atan(sqrt(xx(s1).*6));
s2 = xx>1/6 & xx<=5/6;
f2(s2) = acot2(3/2-xx(s2).*3);
s3 = xx>5/6;
f2(s3) = pi-atan(sqrt(6.*(1-xx(s3))));

xx2 = xx-1/2;
f3 = zeros(1, 200);
absxx = abs(xx2);
s1  = absxx>1/3;
f3(s1) = pi/2-atan(sqrt((1/2-absxx(s1)).*6));
f3(~s1) = pi/2-acot(absxx(~s1).*3);
f3 = f3.*sign(xx2) + pi/2;