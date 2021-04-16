function v = samplesetinobject(u, objecttype, flag_polar)
% P = samplesetinobject(u, objecttype);

if nargin<3
    flag_polar = false;
end

% N = size(u, 1);
switch objecttype
    case 'cube'
        v = u.*2-1;
    case 'sphere'
        rho = nthroot(u(:,1), 3);
        theta = u(:,2).*(pi*2);
        phi = acos(1-u(:,3).*2);
        if flag_polar
            v = [rho theta phi];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    case 'cylinder'
        rho = sqrt(u(:,1));
        theta = u(:,2).*(pi*2);
        z = u(:,3).*2-1;
        if flag_polar
            v = [rho theta z];
        else
            v = [rho.*sin(theta) rho.*cos(theta) z];
        end
    case 'spcylinder'
        theta = u(:,2).*(pi*2);
        phi = invfphi1_spcyl(u(:,3));
        rho = nthroot(u(:,1), 3).*sqrt(min(abs(tan(phi)), abs(cot(phi))).^2+1);
        if flag_polar
            v = [rho.*sin(phi) theta rho.*cos(phi)];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    otherwise
        error('Not support object type %s!', objecttype);
end

end


function phi = invfphi1_spcyl(t)

t = t.*2-1;
abst = abs(t);
s = abst>2/3;
phi = (pi/2 - atan(sqrt(3-3.*abst)).*s - acot(abst.*(3/2)).*(~s)).*sign(t) + pi/2;

end