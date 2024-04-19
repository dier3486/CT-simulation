function [dataflow, prmflow, status] = reconnode_materialzmapprepare(dataflow, prmflow, status)
% prepare node, Z map prepare
% [dataflow, prmflow, status] = reconnode_materialzmapprepare(dataflow, prmflow, status);

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


% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% prms, 
% HU = prmflow.materialdecomp.HU;
hxy = prmflow.recon.voxelsize;
hz = prmflow.recon.imageincrement;
imagesize = int32(prmflow.recon.imagesize);

% to smooth the Z by neighbour
if isfield(nodeprm, 'SmoothRadius')
    SmoothRadius = nodeprm.SmoothRadius;
else
    SmoothRadius = 4.0;
end
if isfield(nodeprm, 'SmoothSigma')
    SmoothSigma = nodeprm.SmoothSigma;
else
    SmoothSigma = 3.0;
end
Cneigh = int32(ceil([SmoothRadius/hxy  SmoothRadius/hxy  SmoothRadius/hz]));
[X, Y, Z] = meshgrid(-Cneigh(1) : Cneigh(1), -Cneigh(2) : Cneigh(2), -Cneigh(3) : Cneigh(3));
X = X(:); Y = Y(:); Z = Z(:);
R = sqrt((single(X).*hxy).^2+(single(Y).*hxy).^2+(single(Z).*hz).^2);
Sr = (R<= SmoothRadius);
% the weight (depending of distance between the voxels) in smoothing Z
prmflow.materialdecomp.SmoothWradius = exp(-R(Sr).^2./SmoothSigma^2./2);
% use to find out the index of the neighbouring voxels
prmflow.materialdecomp.SmoothIndex = Y(Sr) + X(Sr).*imagesize(2) + Z(Sr).*imagesize(2)*imagesize(1);
% NOTE: the SmoothIndex shall be in type int32, while the single or int16 could lay on overflow!

% smoothing Z paramters
if isfield(nodeprm, 'SmoothHalfcut')
    % the half-cut of the density difference in smoothing Z
    prmflow.materialdecomp.SmoothHalfcut = nodeprm.SmoothHalfcut;
else
    prmflow.materialdecomp.SmoothHalfcut = 50;
end
if isfield(nodeprm, 'SmoothSigma')
    % to smooth that half-cut
    prmflow.materialdecomp.SmoothSigma = nodeprm.SmoothSigma;
else
    prmflow.materialdecomp.SmoothSigma = 0.5;
end

% material Z strategy
% parameters in node configure
if isfield(nodeprm, 'AirCut')
    % minimum cut of the density in calculating material B rate
    prmflow.materialdecomp.AirCut = nodeprm.AirCut;
else
    prmflow.materialdecomp.AirCut = 50;
end
if isfield(nodeprm, 'ZlowerBondary')
    % lower limitation of the Z
    prmflow.materialdecomp.ZlowerBondary = nodeprm.ZlowerBondary;
else
%     prmflow.materialdecomp.ZlowerBondary = 2.0;
    prmflow.materialdecomp.ZlowerBondary = prmflow.materialdecomp.ZA*0.1;
end
if isfield(nodeprm, 'ZlowerSigma')
    % to smooth that lower limitation 
    prmflow.materialdecomp.ZlowerSigma = nodeprm.ZlowerSigma;
else
    prmflow.materialdecomp.ZlowerSigma = 0.02;
end
if isfield(nodeprm, 'ZupperBondary')
    % upper limitation of the Z
    prmflow.materialdecomp.ZupperBondary = nodeprm.ZupperBondary;
else
    prmflow.materialdecomp.ZupperBondary = 40.0;
end
if isfield(nodeprm, 'ZupFix')
    % rate in the up-fixing of Z by density due to the natural limitation of the materials' density.
    prmflow.materialdecomp.ZupFix = nodeprm.ZupFix;
else
    % default is 40%
    prmflow.materialdecomp.ZupFix = 0.4;
end
if isfield(nodeprm, 'ZupSigma')
    % to smooth that limitation
    prmflow.materialdecomp.ZupSigma = nodeprm.ZupSigma;
else
    prmflow.materialdecomp.ZupSigma = 0.1;
end
if isfield(nodeprm, 'Zair')
    % the Z of air. We can not measure the Z of air, this value will be used to fill up the blanks.
    prmflow.materialdecomp.Zair = nodeprm.Zair;
else
    prmflow.materialdecomp.Zair = 7.0;
end
% log(9) to be used
prmflow.materialdecomp.log9 = log(9.0);

% color-map consol
if isfield(nodeprm, 'ImageWindow')
    % window of CT-value
    prmflow.materialdecomp.ImageWindow = nodeprm.ImageWindow;
else
    prmflow.materialdecomp.ImageWindow = [100 5000];
end
% the color palette: Hue, Saturation and Lightness
if isfield(nodeprm, 'HueSample')
    % Hue-map sampling
    prmflow.materialdecomp.colorpalette.HueSample = nodeprm.HueSample;
else
    prmflow.materialdecomp.colorpalette.HueSample = 256;
end
if isfield(nodeprm, 'HueLower')
    % Hue lower boundary
    prmflow.materialdecomp.colorpalette.HueLower = nodeprm.HueLower;
else
    % orange
    prmflow.materialdecomp.colorpalette.HueLower = 1/12;
    % red
    % prmflow.materialdecomp.HueLower = 0;
    % yellow
    % prmflow.materialdecomp.HueLower = 1/6;
end
if isfield(nodeprm, 'HueUpper')
    % Hue upper boundary
    prmflow.materialdecomp.colorpalette.HueUpper = nodeprm.HueUpper;
else
    % indigo (blue ~ purple)
    prmflow.materialdecomp.colorpalette.HueUpper = 3/4;
    % blue
    % prmflow.materialdecomp.colorpalette.HueUpper = 2/3;
    % purple
    % prmflow.materialdecomp.colorpalette.HueUpper = 5/6;
end
if isfield(nodeprm, 'HueMid')
    % Hue middle (for water)
    prmflow.materialdecomp.colorpalette.HueMid = nodeprm.HueMid;
else
    % grean
    prmflow.materialdecomp.colorpalette.HueMid = 1/3;
    % cyan
    % prmflow.materialdecomp.HueMid = 1/2;
    % chartreuse (yellow ~ grean)
    % prmflow.materialdecomp.HueUpper = 1/4;
end
if isfield(nodeprm, 'HueLorenz')
    % Lorenz curve of Hue 
    prmflow.materialdecomp.colorpalette.HueLorenz = nodeprm.HueLorenz;
else
    % in (-1, 1). 0 is linear, positive value is 'S' shape, negative value is 'C' shape.
    prmflow.materialdecomp.colorpalette.HueLorenz = 0.5;
end
if isfield(nodeprm, 'HueZlower')
    % lower boundary of Z in hue
    prmflow.materialdecomp.colorpalette.HueZlower = nodeprm.HueZlower;
    % when Z=HueZlower, it will be presented in color HueLower
else
    prmflow.materialdecomp.colorpalette.HueZlower = 2;
end
if isfield(nodeprm, 'HueZupper')
    % upper boundary of Z in hue
    prmflow.materialdecomp.colorpalette.HueZupper = nodeprm.HueZupper;
    % when Z=HueZupper, it will be presented in color HueUpper
else
    prmflow.materialdecomp.colorpalette.HueZupper = 50;
end
if isfield(nodeprm, 'HueZmid')
    % middle of Z in hue
    prmflow.materialdecomp.colorpalette.HueZmid = nodeprm.HueZmid;
    % when Z=HueZmid, it will be presented in color HueMid
else
    prmflow.materialdecomp.colorpalette.HueZmid = prmflow.materialdecomp.ZA;
    % I know ZA is water
end
if isfield(nodeprm, 'Saturation')
    % the saturation
    prmflow.materialdecomp.colorpalette.Saturation = nodeprm.Saturation;
else
    prmflow.materialdecomp.colorpalette.Saturation = 0.5;
end
if isfield(nodeprm, 'LightMin')
    % the minimum limit of light
    prmflow.materialdecomp.colorpalette.LightMin = nodeprm.LightMin;
else
    prmflow.materialdecomp.colorpalette.LightMin = 0.25;
end
if isfield(nodeprm, 'LightLorenz')
    % lightness Lorenz curve
    prmflow.materialdecomp.colorpalette.LightLorenz = nodeprm.LightLorenz;
else
    % in (-1, 1). 0 is linear, positive value is 'J' shape, negative value is 'C' shape (friendly to low density).
    prmflow.materialdecomp.colorpalette.LightLorenz = 0;
    % Note: lightness is gray level.
end
if isfield(nodeprm, 'LightRefWindow')
    % reference window width of that lightness Lorenz curve
    prmflow.materialdecomp.colorpalette.LightRefWindow = nodeprm.LightRefWindow;
else
    prmflow.materialdecomp.colorpalette.LightRefWindow = 4000;
end

% prepare the colorpalette
if isfield(nodeprm, 'showpalette')
    showpalette = nodeprm.showpalette;
else
    showpalette = false;
end
% color palette
[HueZsample, HueMap, LightLorenz, LightMin, Saturation] = ...
    colorpalette(prmflow.materialdecomp.colorpalette, prmflow.materialdecomp.ImageWindow, showpalette);
% to return
prmflow.materialdecomp.HueZsample = HueZsample;
prmflow.materialdecomp.HueMap = HueMap;
prmflow.materialdecomp.LightLorenz = LightLorenz;
prmflow.materialdecomp.LightMin = LightMin;
prmflow.materialdecomp.Saturation = Saturation;

% output color map type: RGB, HSV or HSL
if isfield(nodeprm, 'ColormapType')
    prmflow.materialdecomp.ColormapType = nodeprm.ColormapType;
else
    prmflow.materialdecomp.ColormapType = 'RGB';
end

% double to single
prmflow.materialdecomp = everything2single(prmflow.materialdecomp, 'double', 'single');
% NOTE: do not touch int32.

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [HueZsample, HueMap, LightLorenz, LightMin, Saturation] = colorpalette(CP, ImageWindow, show_onoff)

% Hue
HueZsample = single(linspace(CP.HueZlower, CP.HueZupper, CP.HueSample));
% HueMap
beta = (CP.HueUpper - CP.HueMid)/(CP.HueMid - CP.HueLower) * (CP.HueZmid-CP.HueZlower)/(CP.HueZupper - CP.HueZmid);
HueLorenz_neg = -CP.HueLorenz;
HueLorenz_pos = (beta-1 - (beta+1)*CP.HueLorenz) / (beta+1 -(beta-1)*CP.HueLorenz);
t = HueZsample - CP.HueZmid;
signt = sign(t);
t = abs(t) ./ ((CP.HueZupper - CP.HueZlower)/2 + ((CP.HueZupper + CP.HueZlower)/2 - CP.HueZmid).*signt);
t = min(max(t, 0), 1);
HueLorenz = (HueLorenz_pos + HueLorenz_neg)/2 + (HueLorenz_pos - HueLorenz_neg)/2.*signt;
HueMap = mylorenzcurve(t, HueLorenz) .* ...
    ((CP.HueUpper + CP.HueLower)/2 - CP.HueMid + ((CP.HueUpper - CP.HueLower)/2).*signt);
HueMap = min(max(HueMap + CP.HueMid, 0), 1);

% Light
windowW = ImageWindow(2) - ImageWindow(1);
LightLorenz = (CP.LightRefWindow - windowW + (CP.LightRefWindow + windowW)*CP.LightLorenz) / ...
    (CP.LightRefWindow + windowW + (CP.LightRefWindow - windowW)*CP.LightLorenz);
LightMin = CP.LightMin;

% Saturation
Saturation = CP.Saturation;

if show_onoff
    bandwidth = 64;
    bandL = linspace(0, 1, bandwidth);
    bandLZ = 1 - mylorenzcurve(bandL, LightLorenz).*(1-LightMin);
    w0 = (1000 - ImageWindow(1))/windowW;

    figure();
    title('aa');

    a1 = subplot(3,3,[1 4]);
    Hueband = zeros(CP.HueSample * bandwidth, 3);
    Hueband(:, 1) = repmat(linspace(CP.HueLower, CP.HueUpper,  CP.HueSample)', bandwidth, 1);
    Hueband(:, 2) = Saturation;
    Hueband(:, 3) = repelem(bandL(:), CP.HueSample, 1);

    Hueband = reshape(hsv2rgb(hsl2hsv(Hueband)), CP.HueSample, bandwidth, 3);
    image(rot90(Hueband,2));
    axis equal;
    axis tight;
    axis off;
    
    a4 = subplot(3,3, [8 9]);
    Zband = zeros(CP.HueSample * bandwidth, 3);
    Zband(:, 1) = repelem(HueMap(:), bandwidth, 1);
    Zband(:, 2) = Saturation;
    Zband(:, 3) = repmat(bandLZ(:), CP.HueSample, 1);
    Zband = reshape(hsv2rgb(hsl2hsv(Zband)), bandwidth, CP.HueSample, 3);
    image(flipud(Zband));
    axis equal;
    axis tight;
    axis off;

    a2 = subplot(3,3,[2 3 5 6]);
    plot(HueZsample, HueMap, 'k');
    hold on
    plusC = squeeze(interp1(HueZsample, permute(Zband, [2 1 3]), CP.HueZmid));
    plusC = interp1(bandL, plusC, w0);
    plot(CP.HueZmid, CP.HueMid, '+', 'MarkerSize', 12, 'lineWidth', 2, 'MarkerEdgeColor', plusC);
    axis([min(HueZsample) max(HueZsample) min(HueMap) max(HueMap)]);
    grid on
    Xticks = unique([ceil(min(HueZsample)/10)*10 : 10: floor(max(HueZsample)/10)*10, min(HueZsample)]);
    set(a2, 'Xtick', Xticks);

    a3 = subplot(3,3,7);
    plot(bandL, 1-bandLZ);
    axis equal;
    axis([0 1 0 1]);
    grid on

    drawnow;
end

end

