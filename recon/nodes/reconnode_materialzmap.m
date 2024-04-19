function [dataflow, prmflow, status] = reconnode_materialzmap(dataflow, prmflow, status)
% recon node, after two-material decomposed reconstruction create the Z-map
% [dataflow, prmflow, status] = reconnode_materialzmap(dataflow, prmflow, status);

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_materialzmapprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end

% recon prm
imagesize = int32(prmflow.recon.imagesize);
Nimage = int32(prmflow.recon.Nimage);
% Z map basic prm
HU = prmflow.materialdecomp.HU;
ZA = prmflow.materialdecomp.ZA;
ZB = prmflow.materialdecomp.ZB;
AirCut = prmflow.materialdecomp.AirCut;
% Z map advanced prm
ZlowerBondary = prmflow.materialdecomp.ZlowerBondary;
ZupperBondary = prmflow.materialdecomp.ZupperBondary;
ZlowerSigma = prmflow.materialdecomp.ZlowerSigma;
log9 = prmflow.materialdecomp.log9;
ZupFix = prmflow.materialdecomp.ZupFix;
ZupSigma = prmflow.materialdecomp.ZupSigma;
Zair = prmflow.materialdecomp.Zair;
% Z smooth prm
SmoothIndex = prmflow.materialdecomp.SmoothIndex;
SmoothHalfcut = prmflow.materialdecomp.SmoothHalfcut;
SmoothSigma = prmflow.materialdecomp.SmoothSigma;
SmoothWradius = prmflow.materialdecomp.SmoothWradius;
% color palette rpm
HueZsample = prmflow.materialdecomp.HueZsample;
HueMap = prmflow.materialdecomp.HueMap;
Saturation = prmflow.materialdecomp.Saturation;
LightLorenz = prmflow.materialdecomp.LightLorenz;
LightMin = prmflow.materialdecomp.LightMin;
ImageWindow = prmflow.materialdecomp.ImageWindow;
ColormapType = prmflow.materialdecomp.ColormapType;

% raw image
if GPUonoff
    A0 = gpuArray(reshape(real(dataflow.image), [], Nimage));
    B0 = gpuArray(reshape(imag(dataflow.image), [], Nimage));
else
    A0 = reshape(real(dataflow.image), [], Nimage);
    B0 = reshape(imag(dataflow.image), [], Nimage);
end

% put them to GPU
if GPUonoff
    [Nimage, imagesize, HU, ZA, ZB] = ...
        putinGPU(Nimage, imagesize, HU, ZA, ZB);
    [AirCut, ZlowerBondary, ZupperBondary, ZlowerSigma, log9, ZupFix, ZupSigma, Zair] = ...
        putinGPU(AirCut, ZlowerBondary, ZupperBondary, ZlowerSigma, log9, ZupFix, ZupSigma, Zair);
    [SmoothIndex, SmoothHalfcut, SmoothSigma, SmoothWradius] = ...
        putinGPU(SmoothIndex, SmoothHalfcut, SmoothSigma, SmoothWradius);
    [HueZsample, HueMap, Saturation, LightLorenz, LightMin, ImageWindow] = ...
        putinGPU(HueZsample, HueMap, Saturation, LightLorenz, LightMin, ImageWindow);
end

% to return the rho map
densityAonUH = prmflow.materialdecomp.densityA / prmflow.materialdecomp.HU;
dataflow.imageRhomap = real(dataflow.image) .* densityAonUH;

% Z map 0, by B rate
Zmap = B0./max(A0, AirCut);
Zmap = (max(Zmap.*(ZB^3-ZA^3) + ZA^3, 0)).^(1/3);

% Z map 1, fixed by density
Z1 = min(max(A0, AirCut).*ZA./HU, ZupperBondary);
Wlow = 1./(1+exp(-((Z1 - ZlowerBondary)./ZA).*log9./ZlowerSigma));
Wrho = 1./(1+exp(((Z1 - Zmap)./max(Z1, ZA) - ZupFix).*log9./ZupSigma));
Z1 = (Zmap.*Wrho + Z1.*(1-Wrho)).*Wlow + Zair.*(1-Wlow);

% Z map 2, smooth
Np = imagesize(2)*imagesize(1);
Npall = Np * Nimage;
index_neigh0 = SmoothIndex + (1:Np);
for ii = 1 : Nimage
%     index_ii = (1: Np) + (ii-1).*Np;
    % I know class(ii) is int32
    index_neigh = index_neigh0 + (ii-1).*Np;
    s = (index_neigh < 1) | (index_neigh > Npall);
    index_neigh(s) = 1;

    W = (~s)./(exp((abs(A0(index_neigh) - A0(:, ii)') - SmoothHalfcut).*SmoothSigma)+1);
    W = W.*SmoothWradius;
    Zmap(:, ii) = sum(W.*Z1(index_neigh), 1) ./ sum(W, 1);
end

% to return the Z map
dataflow.imageZmap = gather(reshape(Zmap, imagesize(2), imagesize(1), Nimage));

% image of collor map
imageColor = zeros(Npall, 3, 'like', Zmap);
ZboundU = max(HueZsample);
ZboundL = min(HueZsample);
imageColor(:, 1) = interp1(HueZsample, HueMap, max(min(Zmap(:), ZboundU), ZboundL));
% I know in CUDA the interp1 need not to call the max(min(..)) to employ the close boundary.
imageColor(:, 2) = Saturation;
imageColor(:, 3) = 1 - mylorenzcurve(max(min((A0(:) - ImageWindow(1))./(ImageWindow(2)-ImageWindow(1)), 1), 0), ...
    LightLorenz).*(1-LightMin);

switch ColormapType
    case 'HSL'
        1;
    case 'HSV'
        imageColor = hsl2hsv(imageColor);
    case 'RGB'
        imageColor = hsl2hsv(imageColor);
        imageColor = hsv2rgb(imageColor);
    otherwise
        warning('Unknown color type %s! The HSL color map will be turned.', ColormapType);
end

% to return the colored image
dataflow.imageColor = gather(reshape(imageColor, imagesize(2), imagesize(1), Nimage, 3));

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
