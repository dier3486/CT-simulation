function [dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status)
% prepare node, BP prepare, set FOV, image size, center (XYZ), tilt ... for the images
% [dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status);

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

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% common prepare
prmflow.recon = commonbpprepare(prmflow.recon, prmflow.protocol, prmflow.system, nodeprm);

% recon method
if isfield(nodeprm, 'method') && ~isempty(nodeprm.method)
    recon_method = nodeprm.method;
else
    % default BP method
    switch lower(prmflow.recon.scan)
        case 'axial'
            recon_method = '2D';
        case 'helical'
            recon_method = '';
        otherwise
            recon_method = '';
    end
end
if ~strncmpi(recon_method, prmflow.recon.scan, length(prmflow.recon.scan))
    prmflow.recon.method = [prmflow.recon.scan recon_method];
else
    prmflow.recon.method = recon_method;
end
% The default recon method of axial is Axial2D

% switch recon method
switch lower(prmflow.recon.method)
    case 'axial2d'
        % image center and image number
        prmflow.recon.Nimage = prmflow.recon.Nslice * prmflow.recon.Nshot;
        [prmflow.recon.imagecenter, prmflow.recon.reconcenter_2DBP] = ...
            imagescenterintilt(prmflow.recon.center, prmflow.recon);
    case 'axial3d'
        prmflow.recon = axial3Dprepare(prmflow.recon, nodeprm);
    case {'helical', 'helical3d', 'helicalpiline'}
        % helical is always 3D
        prmflow.recon = helicalprepare(prmflow.recon, nodeprm);
    otherwise
        % do nothing
        1;
        % no topo
end

% check rebin
if isfield(prmflow, 'rebin')
    if ~isfield(prmflow.rebin, 'issloperebin') || ~prmflow.rebin.issloperebin
        if prmflow.protocol.gantrytilt~=0
            % It is a mistake!
            warning(['The reconstruction no longer support previous Axialrebin when gantry tilting! ' ...
                'Please replace the reconnode Axialrebin by the Sloperebin.']);
        end
    end
end

% to single
prmflow.recon = everything2single(prmflow.recon, 'double', 'single');

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpooldata;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function [imagecenter, reconcenter_2DBP] = imagescenterintilt(reconcenter, recon)
% Cin is recon center; Cout is the rotation center on images
% 'small-step' style for tilt shots 

imagecenter = repmat(-reconcenter(:)', recon.Nimage, 1);
% Y shfit
% Yshift = -(recon.imageincrement*tan(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Yshift = -(recon.imageincrement*sin(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Yshift = fliplr(Yshift);
end
Yshift = repmat(Yshift(:), recon.Nshot, 1);
reconcenter_2DBP = imagecenter;
reconcenter_2DBP(:, 2) = reconcenter_2DBP(:, 2) + Yshift;
% Z shift
% Zshift = (recon.imageincrement*sec(recon.gantrytilt)).*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
Zshift = recon.imageincrement.*(-(recon.Nslice-1)/2 : (recon.Nslice-1)/2);
if recon.couchdirection > 0
    Zshift = fliplr(Zshift);
end
Zshift = Zshift(:) - (0:recon.Nshot-1).*(recon.imageincrement.*recon.Nslice).*recon.couchdirection;
Zshift = Zshift(:) - recon.startcouch;
imagecenter = [imagecenter Zshift];
reconcenter_2DBP = [reconcenter_2DBP Zshift];

end


function recon = axial3Dprepare(recon, nodeprm)
% more prepare works for 3D Axial

% Nimage and images center
recon.Nimage = recon.Nslice * recon.Nshot;
[recon.imagecenter, recon.reconcenter_2DBP] = imagescenterintilt(recon.center, recon);

% recon range
reconD = sqrt(sum((recon.FOV/2+abs(recon.center)).^2))*2;
% Neighb and Nextslice
if isfield(nodeprm, 'Neighb')
    Neighb = nodeprm.Neighb;
else
    Rfov = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2);
    Neighb = floor((recon.Nslice*recon.delta_z/2 - (sqrt(recon.SID^2-Rfov^2) - Rfov)/recon.SID*(recon.Nslice-1) ...
             /2*recon.delta_z)/recon.imageincrement) + 2;
end
Neighb = min(Neighb, recon.Nslice/2);
recon.Neighb = Neighb;
recon.Nextslice = recon.Nslice + Neighb*2;

% Zinterp table
if isfield(nodeprm, 'Zinterptablesize')
    tablesize = nodeprm.Zinterptablesize;
else
    tablesize = 512;
end
coneflag = 1;
recon.Zinterp = ZetaEta2TzTable(tablesize, recon.maxFOV, reconD, recon.SID, recon.Nslice, recon.gantrytilt, coneflag);
% I know whether couchdirection the recon.Zinterp is same.

% Z upsampling  matrix or table
recon.Zupsamp = Zupsamplingprepare(recon, nodeprm, 1);

end

function recon = helicalprepare(recon, nodeprm)
% more prepare works for 3D Helical

% imageincrement = recon.imageincrement;
% delta_z = recon.delta_z;
% viewblock is the number of views to loop in each 'block' of raw data
if isfield(nodeprm, 'viewblock') && ~isempty(nodeprm.viewblock)
    viewblock = nodeprm.viewblock;
else
    viewblock = recon.Nviewprot;
end
recon.viewblock = viewblock;

% Z upsampling
recon.Zupsamp = Zupsamplingprepare(recon, nodeprm);

% Cone weight
if isfield(nodeprm, 'ConeWeightScale')
    recon.ConeWeightScale = nodeprm.ConeWeightScale;
else
    recon.ConeWeightScale = 1.0;
end

% imagesnumber per pitch
recon.imagesperpitch = recon.pitchlength/recon.imageincrement;

% the governing of the images related with the views
Reff = min(sqrt(sum(recon.center.^2)) + recon.effFOV/2, recon.maxFOV/2)/recon.SID;
if isfield(nodeprm, 'Nimage')
    Nimage = nodeprm.Nimage;
else
    Nimage = [];
end
[recon.Nimage, recon.Nviewskip, recon.Nviewblock, recon.startimgbyblk, recon.endimgbyblk, recon.startimgbyblk_pi, ...
    recon.endimgbyblk_pi, recon.Zgrid, recon.Zviewshift] = ...
    helicalimagesgovern(Reff, recon.Nview, recon.Nviewprot, recon.Nslice, recon.pitch, recon.imagesperpitch, ...
    recon.viewblock, Nimage);

% image center 
Zshift = (0:recon.Nimage-1).*recon.imageincrement + recon.Nviewskip*recon.pitchlength/recon.Nviewprot ...
         - recon.Zviewshift*recon.imageincrement;
Zshift = -Zshift.*recon.couchdirection - recon.startcouch;
recon.imagecenter = [repmat(-recon.center(:)', recon.Nimage, 1)  Zshift(:)];

end
