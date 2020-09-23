
n = size(A, 1);
m = size(B, 1);
Agpu = gpuArray(A);
Bgpu = gpuArray(B);
viewsgpu = gpuArray(views);
Ztiltgpu = gpuArray(Ztilt);
couchgpu = gpuArray(couch);
Nview = size(viewsgpu, 1);
Oobj = gpuArray(object.O);
invVobj = gpuArray(object.invV);

D = zeros(m, Nview);

% M = zeros(3, 3, Nview, 'gpuArray');
% tic;
% for iview = 1:Nview
%     M(:, :, iview) = rotandtilt(viewsgpu(iview), Ztiltgpu(iview));
% end
% toc;

tic;

% maxview = 1152/16;
maxview = 128;
Nlimit = ceil(Nview/maxview);


for ilim = 1:Nlimit
    v1 = (ilim-1)*maxview + 1;
    if ilim<Nlimit
        v2 = ilim*maxview;
    else
        v2 = Nview;
    end
    Nv = v2-v1+1;
    
    MV = rotandtilt(viewsgpu(v1:v2), Ztiltgpu(v1:v2));
    MV = reshape(permute(MV, [1 3 2]), [], 3) * invVobj;
    MV = reshape(permute(reshape(MV, 3, Nv, 3), [1 3 2]), 3, []);
    OV = (Oobj+couchgpu(v1:v2, :))*invVobj;
    OV = reshape(OV', 1, 3*Nv);
    
    Av = reshape(Agpu*MV - OV, 1, 3, Nv);
    Bv = reshape(Bgpu*MV - OV, m, 3, Nv);
    
%     L = zeros(m, 2, maxview, 'gpuArray');
%     R = zeros(m, 2, maxview, 'gpuArray');
%     % sphere
%     % Lab^2
% %     Lab2 = sum((Av - Bv).^2, 2);
%     L(:, 2, 1:Nv) = sum((Av - Bv).^2, 2);
%     % AxB
%     R(:, 2, 1:Nv) = (Av(:,2,:).*Bv(:,3,:)-Av(:,3,:).*Bv(:,2,:)).^2 + ...
%          (Av(:,3,:).*Bv(:,1,:)-Av(:,1,:).*Bv(:,3,:)).^2 + ...
%          (Av(:,1,:).*Bv(:,2,:)-Av(:,2,:).*Bv(:,1,:)).^2;
%     % d0^2
% %     d02 = S./Lab2;
%     R(:, 1, 1:Nv) = R(:, 2, 1:Nv)./L(:, 2, 1:Nv);
%     % d1;
% %     d1 = sqrt(sum(Av.^2, 2) - d02);
%     L(:, 1, 1:Nv) = sqrt(sum(Av.^2, 2) - R(:, 1, 1:Nv));
%     % d2
% %     d2 = 1 - d02;
% %     d2(d2<0) = 0;
% %     d2 = sqrt(d2);
%     R(:, 2, 1:Nv) = sqrt((1 - R(:, 1, 1:Nv)).*(R(:, 1, 1:Nv)<1));
% %     R(:, 2, 1:Nv) = 1 - R(:, 1, 1:Nv);
% %     R(:, 2, 1:Nv) = R(:, 2, 1:Nv).*(R(:, 2, 1:Nv)>0);
% %     R(:, 2, 1:Nv) =sqrt(R(:, 2, 1:Nv));
%     % Lab
%     L(:, 2, 1:Nv) = sqrt(L(:, 2, 1:Nv));
%     
% %      L = [(d1 - d2)./sqrt(Lab2), zeros(m, 1, Nv)];
% %      R = [(d1 + d2)./sqrt(Lab2), ones(m, 1, Nv)];
%     
%     % R1
%     R(:, 1, 1:Nv) = (L(:, 1, 1:Nv)+R(:, 2, 1:Nv))./L(:, 2, 1:Nv);
%     % L1
%     L(:, 1, 1:Nv) = (L(:, 1, 1:Nv)-R(:, 2, 1:Nv))./L(:, 2, 1:Nv);
%     
%     % R2
%     R(:, 2, 1:Nv) = 1;
%     % L2
%     L(:, 2, 1:Nv) = 0;
%     
%     
%     % done
% %     D(:, v1:v2) = squeeze(gather(multiinsect(L(:,:,1:Nv), R(:,:,1:Nv))));
%     Di = min(R(:,:,1:Nv), [], 2) - max(L(:,:,1:Nv), [], 2);
%     Di = Di.*(Di>0);
%     D(:, v1:v2) = squeeze(gather(Di));


	Di = linesinobjectGPU(Av, Bv, 'cube', [], 1);
%     Di = linesinobject(Av, Bv, 'cube');
    D(:, v1:v2) = squeeze(gather(Di));
%     %2
%     Av = permute(Av, [1 3 2]);
%     Bv = permute(Bv, [1 3 2]);
%     Lab2 = sum((Av - Bv).^2, 3);
%     % AxB
%     S  = (Av(:,:,2).*Bv(:,:,3)-Av(:,:,3).*Bv(:,:,2)).^2 + ...
%          (Av(:,:,3).*Bv(:,:,1)-Av(:,:,1).*Bv(:,:,3)).^2 + ...
%          (Av(:,:,1).*Bv(:,:,2)-Av(:,:,2).*Bv(:,:,1)).^2;
     
end

toc;

% n = size(A, 1);
% m = size(B, 1);
% L = sqrt(sum((repelem(A, m, 1) - repmat(B, n, 1)).^2, 2));
% D = zeros(m*n, Nv, class(L));
% for iview = 1:Nv
%     %             vi = views(iview);
%     %             Mi = [cos(vi)  sin(vi)  0;
%     %                   -sin(vi)  cos(vi)   0;
%     %                   0        0         1];
%     Mi = rotandtilt(views(iview), Ztilt(iview));
%     Oi =  object.O + couch(iview, :);
%     Av = (A*Mi - repmat(Oi, n, 1)) * object.invV;
%     Bv = (B*Mi - repmat(Oi, m, 1)) * object.invV;
%     Av = repelem(Av, m, 1);
%     Bv = repmat(Bv, n, 1);
%     D(:, iview) = linesinobject(Av, Bv, inkey, object.Cimage) .* L;
% end


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