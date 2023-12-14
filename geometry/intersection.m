function [D, L] = intersection(A, B, object, type, views, couch, Ztilt, GPUonoff)
% intersections of lines with objects
% [D, L] = intersection(A, B, object, type, views, couch, Ztilt, GPUonoff)

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
if nargin<7
    Ztilt = 0;
end
if nargin<8
    GPUonoff = false;
end
% in GPU mode the A, B, views, couch, Ztilt shall in GPU array

% views are the view angles
views = views(:);
Nview = size(views, 1);
% couch movements
if size(couch, 1) == 1
    couch = repmat(couch, Nview, 1);
end
if size(couch, 2) == 1
    couch = [zeros(Nview, 2) couch];
end
% tilt
if size(Ztilt, 1) == 1
    Ztilt = repmat(Ztilt, Nview, 1);
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
    case {'cube', 'cuboid', 'parallelepiped', 'parallel hexahedron'}
        inkey = 'cube';
    case 'image2d'
        inkey = 'image2D';
%         isimage = true;
    case {'image3d', 'images'}
        inkey = 'image3D';
%         isimage = true;
    otherwise
        inkey = 'nothing';
end

if GPUonoff
    Aclass = classUnderlying(A);
    object_V = gpuArray(cast(object.vector, Aclass));
    object_invV = gpuArray(cast(object.invV, Aclass));
    object_O = gpuArray(cast(object.O, Aclass));
    crossplane = gpuArray(cast(object.crossplane, Aclass));
else
    Aclass = class(A);
    object_V = object.vector;
    object_invV = object.invV;
    object_O = object.O;
    crossplane = object.crossplane;
end

% cross plane(s)
if ~isempty(crossplane)
    % the planes are coveriate
    crossplane(:, 1:3) = crossplane(:, 1:3) * object_V';
end

switch type
    case 'lines'
        % 'lines' means size(A,1) == size(B,1) or size(A,1) == 1
        L = sqrt(sum((A-B).^2, 2));
        Av = (A - object_O) * object_invV;
        Bv = (B - object_O) * object_invV;
        D = linesinobject(Av, Bv, crossplane, inkey, object.Cimage) .* L;
        D = gather(D);
    case 'net'
        % 'net' means size(A,1)=n, size(B,1)=m and the return D is in (n, m)
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        Av = (A - repmat(object_O, n, 1)) * object_invV;
        Bv = (B - repmat(object_O, m, 1)) * object_invV;
        Av = repelem(Av, m, 1);
        Bv = repmat(Bv, n, 1);
        D = linesinobject(Av, Bv, crossplane, inkey, object.Cimage) .* L;
        D = reshape(gather(D), m, n);
    case {'views', 'views-lines'}
        % 'views' means the type 'lines' in rotation (with couch movment and/or gantry tilt)
        maxview = 128;
        Nlimit = ceil(Nview/maxview);
        n = max(size(A, 1), size(B, 1));
        L = sqrt(sum((A-B).^2, 2));
        D = zeros(n, Nview, Aclass);
        for ilim = 1:Nlimit
            v1 = (ilim-1)*maxview + 1;
            if ilim<Nlimit
                v2 = ilim*maxview;
            else
                v2 = Nview;
            end
            Nv = v2-v1+1;
            % rotation matrix(s)
            MV = rotandtilt(views(v1:v2), Ztilt(v1:v2));
            MV = reshape(permute(MV, [1 3 2]), [], 3) * object_invV;
            MV = reshape(permute(reshape(MV, 3, Nv, 3), [1 3 2]), 3, []);
            % O
            OV = (object_O+couch(v1:v2, :))*object_invV;
            OV = reshape(OV', 1, 3*Nv);
            % Av Bv
            Av = reshape(A*MV - OV, [], 3, Nv);
            Bv = reshape(B*MV - OV, [], 3, Nv);
            % Di
            Di = linesinobject(Av, Bv, crossplane, inkey).*L;
            % to D
            D(:, v1:v2) = squeeze(gather(Di));
        end
    case {'views-net'}
        % 'views-net' means the type 'net' in rotation (with couch movment and/or gantry tilt)
        maxview = 32;
        Nlimit = ceil(Nview/maxview);
        n = size(A, 1);
        m = size(B, 1);
        L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
        D = zeros(m*n, Nview, Aclass);
        for ilim = 1:Nlimit
            v1 = (ilim-1)*maxview + 1;
            if ilim<Nlimit
                v2 = ilim*maxview;
            else
                v2 = Nview;
            end
            Nv = v2-v1+1;
            % rotation matrix(s)
            MV = rotandtilt(views(v1:v2), Ztilt(v1:v2));
            MV = reshape(permute(MV, [1 3 2]), [], 3) * object_invV;
            MV = reshape(permute(reshape(MV, 3, Nv, 3), [1 3 2]), 3, []);
            % O
            OV = (Oobj+couch(v1:v2, :))*object_invV;
            OV = reshape(OV', 1, 3*Nv);
            % Av Bv
            Av = reshape(A*MV - OV, [], 3, Nv);
            Bv = reshape(B*MV - OV, [], 3, Nv);
            Av = repelem(Av, m, 1);
            Bv = repmat(Bv, n, 1);
            % Di
            Di = linesinobject(Av, Bv, crossplane, inkey).*L;
            % to D
            D(:, v1:v2) = reshape(gather(Di), m*n, Nv);
        end
    otherwise
        D = [];
        L = [];
end

end

function M = rotandtilt(viewangle, tiltangle)

% I know
% V = [cos(viewangle)  sin(viewangle)   0;
%     -sin(viewangle)  cos(viewangle)   0;
%      0               0                1];
%  
% T = [1   0                0;
%      0   cos(tiltangle)   sin(tiltangle);
%      0  -sin(tiltangle)   cos(tiltangle)];
% M = V*T;

viewangle = reshape(viewangle, 1, 1, []);
tiltangle = reshape(tiltangle, 1, 1, []);
n = size(viewangle, 3);

M = [cos(viewangle)    sin(viewangle).*cos(tiltangle)   sin(viewangle).*sin(tiltangle);
    -sin(viewangle)    cos(viewangle).*cos(tiltangle)   cos(viewangle).*sin(tiltangle);
     zeros(1,1,n)     -sin(tiltangle)                   cos(tiltangle)                ];

end
 