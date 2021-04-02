function v = samplesetinobject(u, objecttype)
% P = samplesetinobject(u, objecttype);

% N = size(u, 1);
switch objecttype
    case 'cube'
        v = u.*2-1;
    case 'sphere'
        rho = nthroot(u(:,1), 3);
        theta = u(:,2).*(pi*2);
        phi = acos(1-u(:,3).*2);
        v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
    case 'cylinder'
        rho = sqrt(u(:,1));
        theta = u(:,2).*(pi*2);
        z = u(:,3).*2-1;
        v = [rho.*sin(theta) rho.*cos(theta) z];
    otherwise
        error('Not support object type %s!', objecttype);
end

end
