gpuD = gpuDevice;

clear;
Nview = 1152;
N = 512;
Nc = 16;

pb.Np = 900;
pb.delta_d = 0.9;
pb.midchannel = 505.5;
pb.viewangle = linspace(0, 2*pi*(Nview-1)./Nview, Nview);
pb.h = 1.0;

Cimage = rand(N, N, Nc, 'single');
% Cimage = repmat(rand(N, N, 1, 'single'), 1, 1, Nc);
Cclass = class(Cimage);

% % CPU
% tic;
% D1 = parallelprojinimage(pb, Cimage, '2D linearinterp');
% toc;

% GPU
% prepare
CimageGPU = gpuArray([reshape(Cimage, [], Nc); zeros(1, Nc)]);
d_h = gpuArray(cast(((1:pb.Np)'-pb.midchannel).*pb.delta_d./pb.h, Cclass));
D2GPU = zeros(pb.Np*Nc, Nview, Cclass, 'gpuArray');
% pbGPU.Np = gpuArray(pb.Np);
% pbGPU.delta_d = gpuArray(pb.delta_d);
% pbGPU.midchannel = gpuArray(pb.midchannel);
% pbGPU.viewangle = gpuArray(pb.viewangle);
% pbGPU.h = gpuArray(pb.h);


Nx = gpuArray(N);
Ny = gpuArray(N);
Nz = gpuArray(Nc);
iNx = gpuArray(int64(N));
iNy = gpuArray(int64(N));
Np = gpuArray(pb.Np);
Nall = Nx*Ny;
viewangle = gpuArray(cast(pb.viewangle(:), Cclass));
viewindex = gpuArray(1:Nview)';


h = gpuArray(cast(pb.h, Cclass));

phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);
ctgviewangle = tan(pi/2-viewangle);
tanviewangle = tan(viewangle);
cscviewangle = csc(viewangle);
secviewangle = -sec(viewangle);

% % 1
% tic
% for iview = viewindex'
%     [inter_a, index_1, index_2, cs_view] = ...
%         parallellinearinterp2D(Nx, Ny, d_h, viewangle(iview), Nall);
%     D2GPU(:, iview) = reshape(sum(reshape(CimageGPU(index_1(:), :).*(1-inter_a(:)) + ...
%         CimageGPU(index_2(:), :).*inter_a(:), Np, [], Nz), 2).*(abs(cs_view)*h), [], 1);
% end
% D2 = gather(D2GPU);
% toc;

% 2
maxview = 10;
maxview = gpuArray(maxview);
ty = cast(-(Ny-1)/2:(Ny-1)/2, Cclass);
tx = cast(-(Nx-1)/2:(Nx-1)/2, Cclass);
index_t = zeros(Np*maxview, N, Cclass, 'gpuArray');
tic;

% phase1
Nview_phase = sum(phase_flag);
phase_index = find(phase_flag);
Nlimit = ceil(Nview_phase/maxview);
% for iview = viewindex(phase_flag)
for i_lim = 1:Nlimit
    % view number for each step
    if i_lim < Nlimit
        Nview_lim = maxview;
    else
        Nview_lim = Nview_phase - maxview*(Nlimit-1);
    end
    % index of view angles 
    index_lim = phase_index((1:Nview_lim) + maxview*(i_lim-1));

        dt = ctgviewangle(index_lim);
        cs_view = cscviewangle(index_lim);
%         Nx_x = Nx;
%         Ny_y = Ny;

    t = repelem(dt.*ty, Np, 1) - repmat(reshape(d_h.*cs_view', [], 1), 1, Ny) + 1/2;
    index_t = floor(t);
    inter_alpha = t - index_t;
    
    index_1 = index_t + Nx/2;
    index_2 = index_1 + 1;
    s1 = index_1<=0 | index_1>Nx;
    s2 = index_2<=0 | index_2>Nx;
%     index_1(index_1<=0 | index_1>Nx_x) = nan;
%     index_2(index_2<=0 | index_2>Nx_x) = nan;
    
%         index_1 = index_1 + repmat((0:Ny-1).*Nx, Np, 1);
%         index_2 = index_2 + repmat((0:Ny-1).*Nx, Np, 1);
        index_1 = index_1 + (0:Ny-1).*Nx;
        index_2 = index_2 + (0:Ny-1).*Nx;

%     index_1(s1) = Nc+1;
%     index_2(s2) = Nc+1;
    index_1 = index_1.*(~s1) + s1.*(Nall+1);
    index_2 = index_2.*(~s2) + s2.*(Nall+1);
    
%     D2GPU(:, index_lim) = reshape(sum(CimageGPU(index_1).*(1-inter_alpha) + ...
%         CimageGPU(index_2).*inter_alpha, 2), Np, Nview_lim).*(abs(cs_view').*h);
    
    D2GPU(:, index_lim) = reshape(permute(reshape(sum(reshape(CimageGPU(index_1(:), :).*(1-inter_alpha(:)) + ...
        CimageGPU(index_2(:), :).*inter_alpha(:), Np*Nview_lim, Ny, Nz), 2), Np, Nview_lim, Nz), ...
        [1 3 2]), Np*Nz, Nview_lim).*(abs(cs_view')*h);
%     toc;
end

% phase2
Nview_phase = sum(~phase_flag);
phase_index = find(~phase_flag);
Nlimit = ceil(Nview_phase/maxview);
for i_lim = 1:Nlimit
    % view number for each step
    if i_lim < Nlimit
        Nview_lim = maxview;
    else
        Nview_lim = Nview_phase - maxview*(Nlimit-1);
    end
    % index of view angles 
    index_lim = phase_index((1:Nview_lim) + maxview*(i_lim-1));

        dt = tanviewangle(index_lim);
        cs_view = secviewangle(index_lim);
        
%         Nx_x = Ny;
%         Ny_y = Nx;

%     t = repmat((-(Ny_y-1)/2:(Ny_y-1)/2).*dt, Np, 1) - repmat(d_h.*cs_view, 1, Ny_y) + 1/2;
    t = repelem(dt.*tx, Np, 1) - repmat(reshape(d_h.*cs_view', [], 1), 1, Nx) + 1/2;
    
    index_t = floor(t);
    inter_alpha = t - index_t;
    
    index_1 = index_t + Ny/2;
    index_2 = index_1 + 1;
    s1 = index_1<=0 | index_1>Ny;
    s2 = index_2<=0 | index_2>Ny;
%     index_1(index_1<=0 | index_1>Nx_x) = nan;
%     index_2(index_2<=0 | index_2>Nx_x) = nan;
    

        index_1 = (index_1-1).*Nx + (1:Nx);
        index_2 = (index_2-1).*Nx + (1:Nx);

%     index_1(s1) = Nall+1;
%     index_2(s2) = Nall+1;
    index_1 = index_1.*(~s1) + s1.*(Nall+1);
    index_2 = index_2.*(~s2) + s2.*(Nall+1);
    
%     D2GPU(:, iview) = sum(CimageGPU(index_1).*(1-inter_alpha) + ...
%         CimageGPU(index_2).*inter_alpha, 2).*(abs(cs_view)*h);
%     D2GPU(:, index_lim) = reshape(sum(CimageGPU(index_1).*(1-inter_alpha) + ...
%         CimageGPU(index_2).*inter_alpha, 2), Np, Nview_lim).*(abs(cs_view').*h);
    D2GPU(:, index_lim) = reshape(permute(reshape(sum(reshape(CimageGPU(index_1(:), :).*(1-inter_alpha(:)) + ...
        CimageGPU(index_2(:), :).*inter_alpha(:), Np*Nview_lim, Nx, Nz), 2), Np, Nview_lim, Nz), ...
        [1 3 2]), Np*Nz, Nview_lim).*(abs(cs_view')*h);
%     toc;
end
D3 = gather(D2GPU);
toc;



