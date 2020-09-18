gpuD = gpuDevice;

clear;
N = 1000000;
A = rand(N,3).*10;
B = rand(N,3).*10;

objecttype = 'sphere';

tic;
D = linesinobject(A, B, objecttype);
toc;

tic;
Agpu = gpuArray(A);
Bgpu = gpuArray(B);
toc;

tic;
Dgpu = linesinobjectGPU(Agpu, Bgpu, objecttype);
toc;


tic;
LR = zeros(N, 4, 'gpuArray');
% Lab^2
LR(:, 1) = sum((Agpu - Bgpu).^2, 2);
% AxB
LR(:, 3) = (Agpu(:,2).*Bgpu(:,3)-Agpu(:,3).*Bgpu(:,2)).^2 + ...
         (Agpu(:,3).*Bgpu(:,1)-Agpu(:,1).*Bgpu(:,3)).^2 + ...
         (Agpu(:,1).*Bgpu(:,2)-Agpu(:,2).*Bgpu(:,1)).^2;
toc
% d0^2 = AxB / Lab^2
LR(:, 3) = LR(:, 3)./LR(:, 1);
% d1 = sqrt(A^2-d0^2)
LR(:, 2) = sqrt(sum(Agpu.^2, 2) - LR(:, 3));
% d2 = sqrt(1-d0^2);
LR(:, 4) = 1.0 - LR(:, 3);
LR((LR(:, 4)<0), 4) = 0;
LR(:, 4) = sqrt(LR(:, 4));
% Lab
LR(:, 1) = sqrt(LR(:, 1));
LR(:, 3) = LR(:, 1);
% L, R
LR(:, 1) = (LR(:, 2) - LR(:, 4))./LR(:, 1);
LR(:, 3) = (LR(:, 2) + LR(:, 4))./LR(:, 3);
% 0 , 1
LR(:, 2) = 0;
LR(:, 4) = 1;
% D
D2 = multiinsect(LR(:,[1 2]), LR(:, [3 4]));

toc;