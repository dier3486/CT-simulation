function rebin = helicalrebinprepare(rebin, detector, fanangles, focalangle, Nviewprot)
% rebin prepare for Helical
% support XDFS ,not support Z-DFS, no QDO in helical
% recon = rebinprepare(rebin, detector, focalposition, Nview, isQDO);
% where the inputs,
%   detector,                   the struct of detector corr, e.g. prmflow.system.detector,
%   fanangles,                  they are,
%   focalangle,                 [fanangles, focalangle] = detpos2fanangles(detposition, focalposition);                              
%   Nviewprot,                  the view number per rotation, e.g. prmflow.recon.Nviewprot,
% The returns are,
%   rebin.delta_view,           delta view angle
%   rebin.faninterpkern,        interp coeeficients from fan-angles to 'ideal' fan-angles for fan-Radial rebin
%   rebin.dfan,                 equal fan size
%   rebin.idealphi,             equal radial fan-ganles (ideal phi)
%   rebin.Npixel,               Npixel
%   rebin.Nviewprot,            Nviewprot after rebin
%   rebin.Nreb,                 pixel number after rebin
%   rebin.delta_d,              pixle size after rebin
%   rebin.midchannel,           midchannel after rebin
%   rebin.midU,                 midchannel before slope fan-Radial rebin
%   rebin.DFSviewinterp,        DFS

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

[Npixel, Nslice, ~] = size(fanangles);
Nps = Npixel*Nslice;
midU = detector.mid_U;
hx_ISO = detector.hx_ISO;
SID = detector.SID;
Nfocal = length(focalangle);
Npixelmf = Npixel*Nfocal;
% fanangles-focalangle
% fanangles = reshape(reshape(fanangles, Nps, []) - focalangle, Npixel, Nslice, Nfocal);
fanangles = reshape(fanangles, Nps, []) - focalangle;
fanangles = reshape(fanangles', Npixelmf, Nslice);

% mid_U
if Nfocal == 2
    midU = DFSmidchannel(midU(1), abs(focalangle(1)-pi/2) > abs(focalangle(2)-pi/2));
else
    midU = midU(1);
end
% delta d
delta_d = hx_ISO/Nfocal;

% ideal equal fan angle and ideal equal radial phi
[~, dfan, idealphi, midU_phi] = idealfanangles(fanangles, midU, delta_d/SID);
Nreb = length(idealphi);

faninterpkern = zeros(Nreb, Nslice, 'single');
for islice = 1:Nslice
    faninterpkern(:, islice) = interp1(fanangles(:,islice), 1:Npixelmf, idealphi, 'linear', 'extrap');
end

% delta view
delta_view = single(pi*2/Nviewprot);

% Nextraview
Nextraview = [-floor(min(idealphi)/delta_view/Nfocal) ceil(max(idealphi)/delta_view/Nfocal)];
% Only for Azi after Radial.
% We should employ a swicth to distinct Azi after Radial or Radial after Azi, not here.

% DFS
if Nfocal == 2
    DFSviewinterp = -((focalangle-pi/2)./delta_view + [-1/2 1/2])./2;
    DFSviewshift = sum(DFSviewinterp) + 0.5;
else
    DFSviewinterp = [];
    DFSviewshift = 0;
end

% prepare for fan-Radial
rebin.Nviewprot = Nviewprot;    % I know it is raw.Nviewprot.
rebin.Npixel = Npixel;
rebin.Nslice = Nslice;
rebin.Nreb = Nreb;
rebin.Nfocal = Nfocal;
rebin.delta_d = delta_d;
rebin.midchannel = midU;
rebin.midU_phi = midU_phi;
rebin.faninterpkern = faninterpkern;
rebin.dfan = dfan;
rebin.idealphi = idealphi;
rebin.DFSviewinterp = DFSviewinterp;
rebin.DFSviewshift = DFSviewshift;
rebin.delta_view = delta_view;
rebin.gantrytilt = 0;
rebin.Nextraview = Nextraview;

% other prms from detector
rebin.SID = detector.SID;
rebin.delta_z = detector.hz_ISO;


end


