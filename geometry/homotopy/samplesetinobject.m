function v = samplesetinobject(u, objecttype, flag_outpolar, flag_inpolar)
% P = samplesetinobject(u, objecttype);

if nargin<3
    flag_outpolar = false;
end
if nargin<4
    flag_inpolar = false;
end

% N = size(u, 1);
switch objecttype
    case 'cube'
        % basic cubic
        v = u.*2-1;
    case 'sphere'
        % basic sphere
        rho = nthroot(u(:,1), 3);
        theta = u(:,2).*(pi*2);
        phi = acos(1-u(:,3).*2);
        if flag_outpolar
            v = [rho theta phi];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    case 'cylinder'
        % basic cylinder
        rho = sqrt(u(:,1));
        theta = u(:,2).*(pi*2);
        z = u(:,3).*2-1;
        if flag_outpolar
            v = [rho theta z];
        else
            v = [rho.*sin(theta) rho.*cos(theta) z];
        end
    case 'spcylinder'
        % sphere based cylinder, e.q. to fun(fun(u,'sphere'),'sphere2cylinder');
        theta = u(:,2).*(pi*2);
        t = u(:,3).*2-1;
        abst = abs(t);
        s = abst>2/3;
        phi = (pi/2 - atan(sqrt(3-3.*abst)).*s - acot(abst.*(3/2)).*(~s)).*sign(t) + pi/2;
        rho = nthroot(u(:,1), 3).*sqrt((3-3.*abst).*s + (t.*(3/2)).^2.*(~s)+1);
        if flag_outpolar
            v = [rho.*sin(phi) theta rho.*cos(phi)];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    case 'sphere2cylinder'
        if ~flag_inpolar
            u = xyz2polar(u);
        end
        theta = u(:, 2);
        t = cos(u(:, 3));
        abst = abs(t);
        s = abst>2/3;
        phi = (atan(sqrt(3-3.*abst)).*s + acot(abst.*(3/2)).*(~s)-pi/2).*sign(t)+pi/2;
        rho = u(:, 1).*sqrt((3-3.*abst).*s + (t.*(3/2)).^2.*(~s)+1);
        if flag_outpolar
            v = [rho.*sin(phi) theta rho.*cos(phi)];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    case 'cylinder2sphere'
        if ~flag_inpolar
            u = xyz2cylpolar(u);
        end
        tanphi = fillmissing(-u(:,1)./u(:,3), 'constant', 0);
%         tanphi2 = tanphi.^2;
        s = abs(tanphi)<1;
%         phi(s) = acos((1-tanphi(s).^2./3).*sign(u(s,3)));
%         phi(~s) = acos(-(2/3)./tanphi(~s));
        cosphi = fillmissing((1-tanphi.^2./3).*sign(u(:,3)).*s, 'constant', 0) + ...
                 fillmissing(-(2/3)./tanphi.*(~s), 'constant', 0);
        phi = acos(cosphi);
        rho = sqrt(u(:,1).^2+u(:,3).^2);
        theta = u(:,2);
        if flag_outpolar
            v = [rho theta phi];
        else
            v = [rho.*sin(phi).*sin(theta) rho.*sin(phi).*cos(theta) rho.*cos(phi)];
        end
    otherwise
        error('Not support object type %s!', objecttype);
end

end


% function phi = invfphi1_spcyl(t)
% 
% t = t.*2-1;
% abst = abs(t);
% s = abst>2/3;
% phi = (pi/2 - atan(sqrt(3-3.*abst)).*s - acot(abst.*(3/2)).*(~s)).*sign(t) + pi/2;
% 
% end





function v = xyz2cylpolar(u)
% xyz to polar

v = u;
% rho
v(:,1) = sqrt(sum(u(:,1:2).^2, 2));
% theta
v(:,2) = mod(atan2(u(:,2), u(:,1)), pi*2);

end