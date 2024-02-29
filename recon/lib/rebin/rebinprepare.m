function rebin = rebinprepare(rebin, detector, fanangles, focalangle, Nviewprot, isQDO)
% rebin prepare for Axial
% support QDO, XDFS and QDO+XDFS, but not support Z-DFS
% recon = rebinprepare(rebin, detector, focalposition, Nview, isQDO);
% where the inputs,
%   detector,                   the struct of detector corr, e.g. prmflow.system.detector,
%   fanangles,                  they are,
%   focalangle,                [fanangles, focalangle] = detpos2fanangles(detposition, focalposition);                              
%   Nviewprot,                  the view number per rotation, e.g. prmflow.recon.Nviewprot,
%   isQDO,                      the flag (a bool value) of if do QDO;
% The returns are,
%   rebin.delta_view,           delta view angle
%   rebin.interalpha_azi,       interp coeeficients for Azi-rebin
%   rebin.vindex1_azi,
%   rebin.vindex2_azi,          index for Azi-rebin
%   rebin.Npixel,               Npixel or Npixel*2 for QDO
%   rebin.Nviewprot,            Nviewprot or Nviewprot/2 for QDO
%   rebin.QDOorder,             the alternating order of the pixels for QDO
%   rebin.interalpha_rad,       interp coeeficients radial-rebin
%   rebin.radialindex,          index for radial-rebin
%   rebin.Nreb,                 pixel number after rebin
%   rebin.delta_d,              pixle size after rebin
%   rebin.midchannel,           midchannel after rebin

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

% default, not QDO
if nargin<5
    isQDO = false;
end

% parameters to use
% Npixel = double(detector.Npixel);
% Nslice = double(detector.Nmergedslice);
[Npixel, Nslice, ~] = size(fanangles);
% Nslice = 1; % debug
Nps = Npixel*Nslice;
mid_U = detector.mid_U;
hx_ISO = detector.hx_ISO;
SID = detector.SID;
% Nfocal = size(focalposition, 1);
Nfocal = length(focalangle);
% NOTE: focal number is controled by the size of focalposition

% fanangles-focalangle
fanangles = reshape(fanangles, Nps, []) - focalangle;
% I know fanangles was in size Nps*Nfocal or Nps*1

% perpare for Azi rebin
Nviewprot_focal = Nviewprot/Nfocal;
delta_view = pi*2/Nviewprot_focal;
% f is phi/delta_view, where phi is the phi-angle of each pixel in 1st view
f = fanangles./delta_view + (0:Nfocal-1)./Nfocal;
% but for DFS the f contains two columns, for the 1st foal and 2nd focal, 
% the 2nd focal will move one view so need to +0.5, that is (0:Nfocal-1)./Nfocal = [0 0.5];
f = f(:);
% viewindex is the view index of the pixels will be rebined to in Azirebin in first view
viewindex = double(floor(f));
rebin.delta_view = delta_view;
rebin.interalpha_azi = f-viewindex;
% rebin.interalpha_azi = repmat(f-viewindex, Nslice, 1);

viewindex = viewindex + 1;  % +1 due to the floor started from 0
% for DFS, the viewindex will have two columns, for the 1st and 2nd focal, though (:)ed 
% and the value the two columns shall almost same except several indiviually points.

% startvindex is the max of view index+1 and mod by Nviewprot_focal (to return)
rebin.startvindex = mod(max(viewindex), Nviewprot_focal) + 1;

viewindex = repmat(viewindex(:), 1, Nviewprot_focal) + repmat(0:Nviewprot_focal-1, Nps*Nfocal, 1);
pixelindex = reshape(1:Nps*Nfocal, Nfocal, [])';
% I know pixelindex = [1 3 5 7 ... 2 4 6 8 ...] for DFS
rebin.vindex1_azi = mod(viewindex-1, Nviewprot_focal).*(Nps*Nfocal) + repmat(pixelindex(:), 1, Nviewprot_focal);
rebin.vindex2_azi = mod(viewindex, Nviewprot_focal).*(Nps*Nfocal) + repmat(pixelindex(:), 1, Nviewprot_focal);

% prepare for radial rebin
% mid_U
if Nfocal == 2
    mid_U = DFSmidchannel(mid_U(1), abs(focalangle(1)-pi/2) > abs(focalangle(2)-pi/2));
else
    mid_U = mid_U(1);
end
% t
if isQDO
    % QDO order
    [a1, a2] = QDOorder(Npixel*Nfocal, mid_U);
    s1 = ~isnan(a1);
    s2 = ~isnan(a2);
    rebin.Npixel = max([a1, a2]);
    rebin.QDOorder = [a1(:), a2(:)];
    % d0 is the distance from ray to ISO
    d0 = SID.*sin(fanangles)';
    d0 = reshape(d0, [], Nslice);
    % QDO d
    d = nan(rebin.Npixel, Nslice);
    d(a1(s1), :) = d0(s1, :);
    d(a2(s2), :) = -d0(s2, :);
    % delta_t and mid_t
    delta_d = hx_ISO/Nfocal/2;
    mid_t = 0.5;
    % Nview
    rebin.Nviewprot = Nviewprot/Nfocal/2;
else
    rebin.Npixel = Npixel*Nfocal;
    d = SID.*sin(fanangles)';
    d = reshape(d, rebin.Npixel, Nslice);
    delta_d = hx_ISO/Nfocal;
    mid_t = mod(mid_U, 1);
    % Nview
    rebin.Nviewprot = Nviewprot/Nfocal;
end
% radial interp grid
t1 = ceil(min(d(:))/delta_d + mid_t);
t2 = floor(max(d(:))/delta_d + mid_t);
Nreb = t2-t1+1;
tt = ((t1:t2)-mid_t)'.*delta_d;
% interp index
fd = d./delta_d + mid_t;
% dindex = floor(fd) - t1 + 2;
dindex = ceil(fd) - t1 + 1;
dindex(dindex<=0) = 1;
dindex(dindex>Nreb) = Nreb+1;
dindex = dindex + (0:Nslice-1).*(Nreb+1);
tindex = nan(Nreb+1, Nslice);
% tindex(dindex) = 1:rebin.Npixel*Nslice;
tindex(dindex) = repmat((1:rebin.Npixel)', 1, Nslice);
tindex = fillmissing(tindex(1:end-1, :), 'previous');
tindex(tindex>=rebin.Npixel) = rebin.Npixel-1;
tindex = tindex + (0:Nslice-1).*rebin.Npixel;

% DFS
if Nfocal == 2
    DFSviewinterp = -((focalangle-pi/2)./delta_view + [-1/2 1/2])./2;
    DFSviewshift = sum(DFSviewinterp) + 0.5;
else
    DFSviewinterp = [];
    DFSviewshift = 0;
end

% got it
rebin.Nviewprot = Nviewprot;
rebin.Npixel = Npixel;
rebin.Nslice = Nslice;
rebin.Nreb = Nreb;
rebin.Nfocal = Nfocal;
rebin.delta_d = delta_d;
rebin.radialindex = tindex;
rebin.interalpha_rad = (tt - d(tindex))./(d(tindex+1)-d(tindex));
rebin.midU_phi = -t1+1+mid_t;   % rebin output to recon
rebin.midchannel = mid_U;
rebin.DFSviewinterp = DFSviewinterp;
rebin.DFSviewshift = DFSviewshift;
rebin.delta_view = delta_view;

% other prms from detector
rebin.SID = detector.SID;
rebin.delta_z = detector.hz_ISO;

end


