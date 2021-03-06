function D = linesinobjectGPU(A, B, objecttype, Cimage)
% the GPU version of linesinobject()
% D = linesinobjectGPU(A, B, objecttype, Cimage)

if nargin<4
    Cimage = [];
end

% is GPU?
GPUonoff = isa(A, 'gpuArray');
if GPUonoff
    Aclass = classUnderlying(A);
else
    Aclass = class(A);
end
[N, Nview] = size(B);
Nview = Nview/3;

switch objecttype
    case 'sphere'
        % a sphere is |r|<=1.
        % GPU buffer (GPU memory is expensive)
        if GPUonoff
            L = zeros(N, 2, Nview, Aclass, 'gpuArray');
            R = zeros(N, 2, Nview, Aclass, 'gpuArray');
        else
            L = zeros(N, 2, Nview, Aclass);
            R = zeros(N, 2, Nview, Aclass);
        end
        % Lab^2
        L(:, 2, :) = sum((A - B).^2, 2);
        % |AxB|^2
        R(:, 2, :) = pararea2(A, B);
        % d0^2
        R(:, 1, :) = R(:, 2, :)./L(:, 2, :);
        % d1;
        L(:, 1, :) = sqrt(sum(A.^2, 2) - R(:, 1, :));
        % d2
        R(:, 2, :) = sqrt((1 - R(:, 1, :)).*(R(:, 1, :)<1));
        % Lab
        L(:, 2, :) = sqrt(L(:, 2, :));
        % R1
        R(:, 1, :) = (L(:, 1, :)+R(:, 2, :))./L(:, 2, :);
        % L1
        L(:, 1, :) = (L(:, 1, :)-R(:, 2, :))./L(:, 2, :);
        % R2
        R(:, 2, :) = 1;
        % L2
        L(:, 2, :) = 0;
        
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);        
        D = D.*(D>0);

    case 'cylinder'
        % a cylinder is |z|<=1 & x^2+y^2<=1.
        if GPUonoff
            L = zeros(N, 3, Nview, Aclass, 'gpuArray');
            R = zeros(N, 3, Nview, Aclass, 'gpuArray');
        else
            L = zeros(N, 3, Nview, Aclass);
            R = zeros(N, 3, Nview, Aclass);
        end
        % Lxy^2
        L(:, 2, :) = sum((A(:,1:2,:) - B(:,1:2,:)).^2, 2);
        % |AxB|_{xy}^2
        R(:, 2, :) = (A(:,1,:).*B(:,2,:) - A(:,2,:).*B(:,1,:)).^2;
        % d0^2
        R(:, 1, :) = R(:, 2, :)./L(:, 2, :);
        % d1;
        L(:, 1, :) = sqrt(sum(A(:,1:2,:).^2, 2) - R(:, 1, :));
        % d2
        R(:, 2, :) = sqrt((1 - R(:, 1, :)).*(R(:, 1, :)<1));
        % Lab
        L(:, 2, :) = sqrt(L(:, 2, :));
        % R1
        R(:, 1, :) = (L(:, 1, :)+R(:, 2, :))./L(:, 2, :);
        % L1
        L(:, 1, :) = (L(:, 1, :)-R(:, 2, :))./L(:, 2, :);
        % R2 L2
        R(:, 3, :) = cast(B(:, 3, :)<A(:, 3, :), Aclass).*2 - 1;
        R(:, 2, :) = (-R(:, 3, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        L(:, 2, :) = (R(:, 3, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));        
        % R3 L3
        R(:, 3, :) = 1;
        L(:, 3, :) = 0;
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
              
    case 'blade'
        % a blade is 0<=z<=1.
        if GPUonoff
            L = zeros(N, 2, Nview, Aclass, 'gpuArray');
            R = zeros(N, 2, Nview, Aclass, 'gpuArray');
        else
            L = zeros(N, 2, Nview, Aclass);
            R = zeros(N, 2, Nview, Aclass);
        end
        % R1 L1
        R(:, 2, :) = cast(B(:, 3, :)<A(:, 3, :), Aclass);
        R(:, 1, :) = (1 - R(:, 2, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        L(:, 1, :) = (R(:, 2, :) - A(:,3,:))./(B(:,3,:) - A(:,3,:));
        % R2 L2
        R(:, 2, :) = 1;
        L(:, 2, :) = 0;
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
        
    case 'cube'
        % a cube is |x|<=1 & |y|<=1 & |z|<=1.
        if GPUonoff
            L = zeros(N, 4, Nview, Aclass, 'gpuArray');
            R = zeros(N, 4, Nview, Aclass, 'gpuArray');
        else
            L = zeros(N, 4, Nview, Aclass);
            R = zeros(N, 4, Nview, Aclass);
        end
        % s
        R(:, 1:3, :) = cast(B<A, Aclass).*2 - 1;
        % L123, R123
        L(:, 1:3, :) = (R(:, 1:3, :) - A)./(B - A);
        R(:, 1:3, :) = (-R(:, 1:3, :) - A)./(B - A);
        % L4 R4
        R(:, 4, :) = 1;
        L(:, 4, :) = 0;
        % multiinsect
        D = min(R, [], 2) - max(L, [], 2);
        D = D.*(D>0);
        
    case 'image2D'
        % 2D image is an image copied on z direction
        % TBC
        error('Not support yet!');
    
    case {'image3D', 'images'}
        % 3D image is an array of images on z direction
        % TBC
        errror('Not support yet!');
        
    otherwise
        D = zeros(N, 1, Aclass);
        return 
end

end