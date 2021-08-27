function D = parallelprojinimageGPU(parallelbeam, Cimage)
% parallel projection on image(s), GPU plus
% D = parallelprojinimage(parallelbeam, Cimage);

% only support method '2D linearinterp'

% ini GPU Device
gpuD = gpuDevice;

% get inputs
if isa(Cimage, 'gpuArray')
    % I know the parallelbeam.* ara all gpuArray
    [Cclass, delta_d, Np, Nx, Ny, Nz, mid_chn, h, viewangle, Nview, maxview] = getinputsGPU(parallelbeam, Cimage);
    Cimage = reshape(Cimage, [], Nz);
else
    % put on GPU
    [Cclass, delta_d, Np, Nx, Ny, Nz, mid_chn, h, viewangle, Nview, maxview] = getinputsCPU(parallelbeam, Cimage);
    Cimage = gpuArray([reshape(Cimage, [], Nz); zeros(1, Nz)]);
end

Nxy = Nx*Ny;
d_h = ((1:Np)' - mid_chn).*delta_d./h;

% ini result
D = zeros(Np*Nz, Nview, Cclass, 'gpuArray');

% tan, ctg, csc, sec and phase of the viewangle
phase_flag = (viewangle>pi/4 & viewangle<=pi*3/4) | (viewangle>pi*5/4 & viewangle<=pi*7/4);
ctgviewangle = tan(pi/2-viewangle);
tanviewangle = tan(viewangle);
cscviewangle = csc(viewangle);
secviewangle = -sec(viewangle);

ty = cast(-(Ny-1)/2:(Ny-1)/2, Cclass);
tx = cast(-(Nx-1)/2:(Nx-1)/2, Cclass);
% index_t = zeros(Np*maxview, Ny, Cclass, 'gpuArray');

% phase1
Nview_phase = sum(phase_flag);
phase_index = find(phase_flag);
Nlimit = ceil(Nview_phase/maxview);
% loop
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
    t = repelem(dt.*ty, Np, 1) - repmat(reshape(d_h.*cs_view', [], 1), 1, Ny) + 1/2;
    index_t = floor(t);
    inter_alpha = t - index_t;
    
    index_1 = index_t + Nx/2;
    index_2 = index_1 + 1;
    s1 = index_1<=0 | index_1>Nx;
    s2 = index_2<=0 | index_2>Nx;

    index_1 = index_1 + (0:Ny-1).*Nx;
    index_2 = index_2 + (0:Ny-1).*Nx;
    index_1 = index_1.*(~s1) + s1.*(Nxy+1);
    index_2 = index_2.*(~s2) + s2.*(Nxy+1);
    
    D(:, index_lim) = reshape(permute(reshape(sum(reshape(Cimage(index_1(:), :).*(1-inter_alpha(:)) + ...
        Cimage(index_2(:), :).*inter_alpha(:), Np*Nview_lim, Ny, Nz), 2), Np, Nview_lim, Nz), ...
        [1 3 2]), Np*Nz, Nview_lim).*(abs(cs_view')*h);
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
    t = repelem(dt.*tx, Np, 1) - repmat(reshape(d_h.*cs_view', [], 1), 1, Nx) + 1/2;
    index_t = floor(t);
    inter_alpha = t - index_t;
    
    index_1 = index_t + Ny/2;
    index_2 = index_1 + 1;
    s1 = index_1<=0 | index_1>Ny;
    s2 = index_2<=0 | index_2>Ny;
   

    index_1 = (index_1-1).*Nx + (1:Nx);
    index_2 = (index_2-1).*Nx + (1:Nx);
    index_1 = index_1.*(~s1) + s1.*(Nxy+1);
    index_2 = index_2.*(~s2) + s2.*(Nxy+1);
    
    D(:, index_lim) = reshape(permute(reshape(sum(reshape(Cimage(index_1(:), :).*(1-inter_alpha(:)) + ...
        Cimage(index_2(:), :).*inter_alpha(:), Np*Nview_lim, Nx, Nz), 2), Np, Nview_lim, Nz), ...
        [1 3 2]), Np*Nz, Nview_lim).*(abs(cs_view')*h);
end

end


function [Cclass, delta_d, Np, Nx, Ny, Nz, mid_chn, h, viewangle, Nview, maxview] = getinputsCPU(parallelbeam, Cimage)
% get the inputs and trans to gpu

% class of the Cimage (single or double)
Cclass = class(Cimage);
% Np is the number of the parallel beams (channel number)
Np = gpuArray(cast(parallelbeam.Npixel,Cclass));
% delta_d is the step length of the parallel beams
delta_d = gpuArray(cast(parallelbeam.delta_d, Cclass));
% Nx Ny Nz is the size of the Cimage
if isfield(parallelbeam, 'Nx')
    Nx = gpuArray(cast(parallelbeam.Nx, Cclass));
    Ny = gpuArray(cast(parallelbeam.Ny, Cclass));
    if isfield(parallelbeam, 'Nz')
        Nz = gpuArray(cast(parallelbeam.Nz, Cclass));
    else
        Csize = size(CimageGPU);
        Nz = gpuArray(cast(Csize(end), Cclass));
    end
else
    [Nx, Ny, Nz] = size(Cimage);
    Nx = gpuArray(cast(Nx, Cclass));
    Ny = gpuArray(cast(Ny, Cclass));
    Nz = gpuArray(cast(Nz, Cclass));
end
% midchannel is the midchennel index of the parallel beams
if isfield(parallelbeam, 'midchannel')
    mid_chn = gpuArray(cast(parallelbeam.midchannel, Cclass));
else
    mid_chn = (Np+1)/2;
end
% h is the voxel size of the image
if isfield(parallelbeam, 'h')
    h = gpuArray(cast(parallelbeam.h, Cclass));
else
    h = ones(1, Cclass, 'gpuArray');
end
% viewangle are the angles of the parallel beams
if isfield(parallelbeam, 'viewangle')
    viewangle =  gpuArray(cast(parallelbeam.viewangle(:), Cclass));
elseif isfield(parallelbeam, 'Nview')
    viewangle = gpuArray(cast(linspace(0, pi*2*(parallelbeam.Nview-1)/parallelbeam.Nview)', Cclass));
else
    viewangle = zeros(1, Cclass, 'gpuArray');
end
% Nview is the number of viewangles
Nview = gpuArray(cast(size(viewangle(:), 1), Cclass));
% maxview is the limit of the view numbers in once loop of projection (to control the memory using)
if isfield(parallelbeam, 'maxview')
    maxview = gpuArray(cast(parallelbeam.maxview, Cclass));
else
    maxview = gpuArray(single(10));
end
    
end

function [Cclass, delta_d, Np, Nx, Ny, Nz, mid_chn, h, viewangle, Nview, maxview] = getinputsGPU(parallelbeam, Cimage)

Cclass = classUnderlying(Cimage);
Np = parallelbeam.Npixel;
delta_d = parallelbeam.delta_d;
Nx = parallelbeam.Nx;
Ny = parallelbeam.Ny;
if isfield(parallelbeam, 'Nz')
    Nz = parallelbeam.Nz;
else
    Csize = size(CimageGPU);
    Nz = gpuArray(Csize(end));
end
if isfield(parallelbeam, 'midchannel')
    mid_chn = parallelbeam.midchannel;
else
    mid_chn = (Np+1)/2;
end
if isfield(parallelbeam, 'h')
    h = parallelbeam.h;
else
    h = ones(1, Cclass, 'gpuArray');
end
if isfield(parallelbeam, 'viewangle')
    viewangle = parallelbeam.viewangle(:);
elseif isfield(parallelbeam, 'Nview')
    viewangle = linspace(0, pi*2*(parallelbeam.Nview-1)/parallelbeam.Nview)';
else
    viewangle = zeros(1, Cclass, 'gpuArray');
end
Nview = gpuArray(cast(size(viewangle(:), 1), Cclass));
if isfield(parallelbeam, 'maxview')
    maxview = parallelbeam.maxview;
else
    maxview = gpuArray(single(10));
end
end