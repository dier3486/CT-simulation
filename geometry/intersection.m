function [D, L] = intersection(A, B, object, type, views, bed)

if nargin<4
    % default type is lines
    type = 'lines';
end
if nargin<5
    views = 0;
end
if nargin<6
    bed = 0;
end

views = views(:);
Nv = size(views, 1);
if size(bed, 1) == 1
    bed = repmat(bed, Nv, 1);
end
if size(bed, 2) == 1
    bed = [zeros(Nv, 2) bed];
end

switch object.type
    case {'sphere', 'ellipsoid'}
        inkey = 'sphere';
    case 'cylinder'
        inkey = 'cylinder';
    case 'blade'
        inkey = 'blade';
    otherwise
        inkey = 'nothing';
end

switch type
    case 'lines'
        n = size(A, 1);
        L = sqrt(sum((A-B).^2, 2));
        Av = (A - repmat(object.O, n, 1)) * object.invV;
        Bv = (B - repmat(object.O, n, 1)) * object.invV;
        D = linesinobject(Av, Bv, inkey) .* L;
    case 'ray'
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        Av = (A - repmat(object.O, n, 1)) * object.invV;
        Bv = (B - repmat(object.O, m, 1)) * object.invV;
        Av = repelem(Av, m, 1);
        Bv = repmat(Bv, n, 1);
        D = linesinobject(Av, Bv, inkey) .* L;
        D = reshape(D, m, n);
    case {'views-lines'}
        n = size(A, 1);
        L = sqrt(sum((A-B).^2, 2));
        D = zeros(n, Nv);
        for iview = 1:Nv
            vi = views(iview);
            Mi = [cos(vi)  -sin(vi)  0;
                  sin(vi)  cos(vi)   0;
                  0        0         1];
            Oi =  object.O + bed(iview, :);
            Av = (A*Mi - repmat(Oi, n, 1)) * object.invV;
            Bv = (B*Mi - repmat(Oi, n, 1)) * object.invV;
            D(:, iview) = linesinobject(Av, Bv, inkey) .* L;
        end
    case {'views', 'views-ray'}
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        D = zeros(m*n, Nv);
        for iview = 1:Nv
            vi = views(iview);
            Mi = [cos(vi)  -sin(vi)  0;
                  sin(vi)  cos(vi)   0;
                  0        0         1];
            Oi =  object.O + bed(iview, :);
            Av = (A*Mi - repmat(Oi, n, 1)) * object.invV;
            Bv = (B*Mi - repmat(Oi, m, 1)) * object.invV;
            Av = repelem(Av, m, 1);
            Bv = repmat(Bv, n, 1);
            D(:, iview) = linesinobject(Av, Bv, inkey) .* L;
        end
        
    otherwise
        D = [];
        L = [];
end

return