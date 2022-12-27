function [D, mu] = projectinphantom(focalposition, detectorposition, phantom, samplekeV, viewangle, couch, gantrytilt, GPUonoff)
% project in phantom(s)
% [D, mu] = projectinphantom(focalposition, detectorposition, phantom, samplekeV, viewangle, couch, gantrytilt, GPUonoff);
% We know Dmu = D*mu; but which could lay out of memory.

% default viewangle(s) and couch
if nargin<5 || isempty(viewangle)
    viewangle = 0;
    couch = 0;
end
if nargin<6 || isempty(couch)
    couch = zeros(size(viewangle));
end
if nargin<7 || isempty(gantrytilt)
    gantrytilt = zeros(size(viewangle));
end
if nargin<8
    GPUonoff = false;
end

if isempty(phantom)
    % do nothing
    D = [];
    mu = [];
    return;
end

% size
Np = size(detectorposition, 1);
Nview = length(viewangle(:));
Nsample = length(samplekeV(:));
Nfocal = size(focalposition, 1);
% viewangle = reshape(viewangle, Nfocal, []);
% couch = reshape(couch, Nfocal, [], 3);
Nobject = phantom.Nobject;

% to GPU
if GPUonoff
    Cclass = 'single';
    focalposition = gpuArray(cast(focalposition, Cclass));
    detectorposition = gpuArray(cast(detectorposition, Cclass));
%     viewangle = reshape(viewangle, Nfocal, []);
    viewangle = gpuArray(cast(viewangle, Cclass));
    couch = gpuArray(cast(couch, Cclass));
    gantrytilt = gpuArray(cast(gantrytilt, Cclass));
else
    Cclass = class(detectorposition);
end
% ini D & mu
D = zeros(Np, Nview, Nobject, Cclass);
mu = zeros(Nobject, Nsample, Cclass);

for iobj = 1:phantom.Nobject
    % loop the objects
    parentobj = phantom.object_tree(iobj);
    object_i = phantom.object{iobj};
    % mu
    mu_i = interp1(object_i.material.samplekeV, object_i.material.mu_total, samplekeV);
    % mu fix
    if parentobj>0
        mu_parent = interp1(phantom.object{parentobj}.material.samplekeV, ...
            phantom.object{parentobj}.material.mu_total, samplekeV);
        mu_i = mu_i - mu_parent;
    end
    mu(iobj, :) = mu_i;
    % L = zeros(Np, Nfocal);
    % I know the L has been done in flewoverbowtie.m
    % fly focal
    for ifocal = 1:Nfocal
        % geometry projection in object
        [D(:, ifocal:Nfocal:end, iobj), ~] = intersection(focalposition(ifocal, :), detectorposition, object_i, 'views', ...
            viewangle(ifocal:Nfocal:end), couch(ifocal:Nfocal:end, :), gantrytilt(ifocal:Nfocal:end), GPUonoff);
    end
end
D = reshape(D, Np, Nview*Nobject);
% I know the Dmu = D*mu; but the outer product could cost too much memory. 

end
