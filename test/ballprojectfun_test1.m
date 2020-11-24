function p = ballprojectfun_test1(viewangle, SID, SDD, Vin)
% test function for ball porjection fitting

[Rb, Zb, phib, x0 , z0, alpha0] = tac(Vin);

x = Rb.*sin(viewangle + phib);
y = Rb.*cos(viewangle + phib);

x = x.*SDD./(SID+y);
z = (SDD*Zb)./(SID+y);

Aalpha = [cos(alpha0) -sin(alpha0); sin(alpha0) cos(alpha0)];
p = Aalpha*[x; z] + [x0; z0];

end