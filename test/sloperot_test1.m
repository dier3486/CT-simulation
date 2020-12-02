Nv = 200;

theta_v = linspace(0, pi*2, Nv);
theta_t = 0.3;

r0 = [0 -500, 0];
r1 = [100, 300, 20];

u0 = zeros(Nv, 3);
u1 = zeros(Nv, 3);
% v0 = zeros(Nv, 3);
% v1 = zeros(Nv, 3);

Rt = [1  0             0           ;
      0  cos(theta_t) -sin(theta_t);
      0  sin(theta_t)  cos(theta_t)];
At = [1  0             0;
      0  sec(theta_t)  0;
      0 -tan(theta_t)  1];
for ii = 1:Nv
    Rv = [cos(theta_v(ii)) -sin(theta_v(ii)) 0;
          sin(theta_v(ii))  cos(theta_v(ii)) 0;
          0                 0                1];

    u0(ii, :) = r0*Rv';
    u1(ii, :) = r1*Rv';
        
end