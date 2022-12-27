function [inter_alpha1, inter_alpha2, index_1, index_2, cs_vangle] = linesinimageLI2D2(Nx, Ny, d_h, theta, AO, Lxy, Next)
% insections of lines in grid-cells image, 2D. linear interp, 1-grid model, finite lines
% [dt, Vindex] = linesinimageLI2D2(Nx, Ny, d_h, theta, AO, L)
% then D =  sum(Cimage(index_1).*inter_alpha1 + Cimage(index_2).*inter_alpha2, 2).*abs(cs_vangle);
% remember to add a 0 after Cimage that Cimage = [Cimage(:); 0]
% WARN: known bug: the Cimage shall be Cimage.'

if nargin<7
    Next = Nx*Ny+1;
end

if isa(d_h, 'gpuArray')
    GPUonoff = true;
    dclass = classUnderlying(d_h);
else
    GPUonoff = false;
    dclass = class(d_h);
end

% mod theta with 2pi
theta = mod(theta(:), pi*2);
% I know theta = theta(:); d = d(:);

% ini
Np = size(theta, 1);
N = max(Nx, Ny);
if GPUonoff
    cs_vangle = gpuArray(zeros(Np, 1, dclass));
    inter_alpha1 = gpuArray(zeros(Np, N, dclass));
    inter_alpha2 = gpuArray(zeros(Np, N, dclass));
    index_1 = gpuArray(ones(Np, N, dclass).*Next);
    index_2 = gpuArray(ones(Np, N, dclass).*Next);
    Nx = gpuArray(Nx);
    Ny = gpuArray(Ny);
else
    cs_vangle = zeros(Np, 1, dclass);
    inter_alpha1 = zeros(Np, N, dclass);
    inter_alpha2 = zeros(Np, N, dclass);
    index_1 = ones(Np, N, dclass).*Next;
    index_2 = ones(Np, N, dclass).*Next;
end

% cos and sin
costheta = cos(theta);
sintheta = sin(theta);

% loop 2 phase
phase_flag = (theta>pi/4 & theta<=pi*3/4) | (theta>pi*5/4 & theta<=pi*7/4);
% phase 1 is the lines close to axis Y; phase 2 is the lines close to X
for iphase = [1 0]
    if iphase
        d_i = d_h(phase_flag);
        dt_i = costheta(phase_flag)./sintheta(phase_flag);
        cs_vangle_i = 1./sintheta(phase_flag);
        L = -AO(phase_flag).*costheta(phase_flag) + d_i.*sintheta(phase_flag);
        R = L + Lxy(phase_flag).*costheta(phase_flag);
        Nx_x = Nx;
        Ny_y = Ny;
        % to return
        cs_vangle(phase_flag) = cs_vangle_i;
    else
        d_i = d_h(~phase_flag);
        dt_i = sintheta(~phase_flag)./costheta(~phase_flag);
        cs_vangle_i = -1./costheta(~phase_flag);
        L = -AO(~phase_flag).*sintheta(~phase_flag) - d_i.*costheta(~phase_flag);
        R = L + Lxy(~phase_flag).*sintheta(~phase_flag);
        Nx_x = Ny;
        Ny_y = Nx;
        % to return
        cs_vangle(~phase_flag) = cs_vangle_i;
    end
    
    t = dt_i*(-(Ny_y-1)/2:(Ny_y-1)/2) - repmat(d_i.*cs_vangle_i, 1, Ny_y) + 1/2;
    index_1_i = floor(t);
    inter_alpha_i = t - index_1_i;
    Lt = t - (L-1/2);
    Lt(Lt<0) = 0;
    Lt(Lt>1) = 1;
    Rt = (R+1/2) - t;
    Rt(Rt<0) = 0;
    Rt(Rt>1) = 1;
    LR = Lt+Rt-1;
%     LR = 1;

    index_1_i = index_1_i + Nx_x/2;
    index_2_i = index_1_i + 1;
    s1 = index_1_i<=0 | index_1_i>Nx_x;
    s2 = index_2_i<=0 | index_2_i>Nx_x;
    if iphase
        index_1_i = index_1_i + (0:Ny-1).*Nx;
        index_2_i = index_2_i + (0:Ny-1).*Nx;
    else
        index_1_i = (index_1_i-1).*Nx + (1:Nx);
        index_2_i = (index_2_i-1).*Nx + (1:Nx);
    end
    index_1_i(s1) = Next;
    index_2_i(s2) = Next;
    
    if iphase
        inter_alpha1(phase_flag, 1:Ny) = (1-inter_alpha_i).*LR;
        inter_alpha2(phase_flag, 1:Ny) = inter_alpha_i.*LR;
        index_1(phase_flag, 1:Ny) = index_1_i;
        index_2(phase_flag, 1:Ny) = index_2_i;
    else
        inter_alpha1(~phase_flag, 1:Nx) = (1-inter_alpha_i).*LR;
        inter_alpha2(~phase_flag, 1:Ny) = inter_alpha_i.*LR;
        index_1(~phase_flag, 1:Nx) = index_1_i;
        index_2(~phase_flag, 1:Nx) = index_2_i;
    end
    
end

return