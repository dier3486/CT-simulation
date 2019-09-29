function [D, L] = intersection(A, B, object, type, views, couch)
% intersections of lines with objects
% [D, L] = intersection(A, B, object, type, views, couch)

% default inputs
if nargin<4
    % default type is lines
    type = 'lines';
end
if nargin<5
    views = 0;
end
if nargin<6
    couch = 0;
end
% views are the view angles
views = views(:);
Nv = size(views, 1);
% couch movements
if size(couch, 1) == 1
    couch = repmat(couch, Nv, 1);
end
if size(couch, 2) == 1
    couch = [zeros(Nv, 2) couch];
end
% % flag for is image
% isimage = false;
% swithc the object type
switch lower(object.type)
    case {'sphere', 'ellipsoid'}
        inkey = 'sphere';
    case 'cylinder'
        inkey = 'cylinder';
    case 'blade'
        inkey = 'blade';
    case 'image2d'
        inkey = 'image2D';
%         isimage = true;
    case {'image3d', 'images'}
        inkey = 'image3D';
%         isimage = true;
    otherwise
        inkey = 'nothing';
end

switch type
    case 'lines'
        n = size(A, 1);
        L = sqrt(sum((A-B).^2, 2));
        Av = (A - repmat(object.O, n, 1)) * object.invV;
        Bv = (B - repmat(object.O, n, 1)) * object.invV;
        D = linesinobject(Av, Bv, inkey, object.Cimage) .* L;
    case 'ray'
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        Av = (A - repmat(object.O, n, 1)) * object.invV;
        Bv = (B - repmat(object.O, m, 1)) * object.invV;
        Av = repelem(Av, m, 1);
        Bv = repmat(Bv, n, 1);
        D = linesinobject(Av, Bv, inkey, object.Cimage) .* L;
        D = reshape(D, m, n);
    case {'views-lines'}
        n = size(A, 1);
        L = sqrt(sum((A-B).^2, 2));
        D = zeros(n, Nv, class(L));
        for iview = 1:Nv
            vi = views(iview);
            Mi = [cos(vi)  sin(vi)  0;
                  -sin(vi)  cos(vi)   0;
                  0        0         1];
            Oi =  object.O + couch(iview, :);
            Av = (A*Mi - repmat(Oi, n, 1)) * object.invV;
            Bv = (B*Mi - repmat(Oi, n, 1)) * object.invV;
            D(:, iview) = linesinobject(Av, Bv, inkey, object.Cimage) .* L;
        end
    case {'views', 'views-ray'}
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        D = zeros(m*n, Nv, class(L));
        for iview = 1:Nv
            vi = views(iview);
            Mi = [cos(vi)  sin(vi)  0;
                  -sin(vi)  cos(vi)   0;
                  0        0         1];
            Oi =  object.O + couch(iview, :);
            Av = (A*Mi - repmat(Oi, n, 1)) * object.invV;
            Bv = (B*Mi - repmat(Oi, m, 1)) * object.invV;
            Av = repelem(Av, m, 1);
            Bv = repmat(Bv, n, 1);
            D(:, iview) = linesinobject(Av, Bv, inkey, object.Cimage) .* L;
        end
        
    otherwise
        D = [];
        L = [];
end

return