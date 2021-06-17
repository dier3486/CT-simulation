function v = polar2xyz(u)
% polar to xyz, u=[rho theta phi], v=[x y z];

v = [u(:,1).*cos(u(:,2)).*sin(u(:,3)) u(:,1).*sin(u(:,2)).*sin(u(:,3)) u(:,1).*cos(u(:,3))];

end