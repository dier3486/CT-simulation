function v = xyz2polar(u)
% xyz to polar, u=[x y z], v=[rho theta phi]. 

v = u;
% rho
v(:,1) = sqrt(sum(u.^2, 2));
% theta
v(:,2) = mod(atan2(u(:,2), u(:,1)), pi*2);
% phi
v(:,3) = fillmissing(acos(u(:,3)./v(:,1)), 'constant', 0);

end