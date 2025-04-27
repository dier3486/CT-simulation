function [dataflow, prmflow, status] = reconnode_helicalbackprojectionrestart(dataflow, prmflow, status)
% helical BP restart
% [dataflow, prmflow, status] = reconnode_helicalbackprojectionrestart(dataflow, prmflow, status);

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

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% update .recon parameters
recon = prmflow.recon;
% I know the recon.Nview was updated in rebin restart

% Nviewact, Zviewshift, Nimage, Zgrid
Nviewact = min(recon.Nview + ceil(recon.viewextra_full), recon.viewnumber - recon.Nviewskip);
Zviewshift = recon.Zviewshift + recon.viewread*recon.imagesperpitch/recon.Nviewprot - recon.availwritten;
Nimage_all = recon.imagewritten - recon.availwritten + round(recon.Nviewact/recon.Nviewprot*recon.imagesperpitch);
imageZgrid = 0 : (Nimage_all-1);

% start/end view of each image in pi-line recon
viewstart_pi = ceil((imageZgrid - Zviewshift).*(recon.Nviewprot/recon.imagesperpitch) - recon.viewextra_pi);
viewend_pi = floor((imageZgrid - Zviewshift).*(recon.Nviewprot/recon.imagesperpitch) + recon.viewextra_pi);

% available images are which can be reconstructed by pi-line
imgavail = viewend_pi <= Nviewact;
Nimage = sum(imgavail);

imageZgrid = imageZgrid(imgavail); 
% re-set the available start/end view
viewstart_pi = viewstart_pi(imgavail);
viewend_pi = viewend_pi(imgavail);

% start/end view of each image in full recon
viewstart_full = ceil((imageZgrid - Zviewshift).*(recon.Nviewprot/recon.imagesperpitch) - recon.viewextra_full);
viewend_full = floor((imageZgrid - Zviewshift).*(recon.Nviewprot/recon.imagesperpitch) + recon.viewextra_full);

% image center 
% Zshift = (0:Nimage-1).*recon.imageincrement + (recon.Nviewskip + recon.viewread)*recon.pitchlength/recon.Nviewprot ...
%          - Zviewshift*recon.imageincrement;
% Zshift = -Zshift.*recon.couchdirection - recon.startcouch;
imageStart = recon.imageStart - recon.availwritten.*recon.couchdirection;
Zshift = (imageStart - (0:Nimage-1).*recon.couchdirection).*recon.imageincrement - recon.startcouch;
prmflow.recon.imagecenter = [repmat(-recon.center(:)', Nimage, 1)  Zshift(:)];
prmflow.recon.imageStart = imageStart;

% InstanceNumber start
prmflow.recon.InstanceStart = recon.InstanceStart + recon.availwritten;

% to return
% Nviewact
prmflow.recon.Nviewact = Nviewact;
% ZviewRange
prmflow.recon.ZviewRange = [0 Nviewact-1];
% start/end view and other
prmflow.recon.viewbyimages_pi = [viewstart_pi; viewend_pi];
prmflow.recon.viewbyimages_full = [viewstart_full; viewend_full];
prmflow.recon.imageZgrid = imageZgrid;
prmflow.recon.Zviewshift = Zviewshift;
% update viewnumber
prmflow.recon.viewnumber = recon.viewnumber - recon.viewread;
prmflow.recon.viewread = 0;

% update Nimage
prmflow.recon.Nimage = Nimage;
prmflow.recon.imagewritten = recon.imagewritten - recon.availwritten;
prmflow.recon.availwritten = 0;

% if pipeline_onoff
%     nodename = status.nodename;
%     nextnode = status.pipeline.(nodename).nextnode;
%     if ~isempty(nextnode) && ~strcmpi(nextnode, 'NULL')
%         newNimage = Nimage - recon.Nimage + recon.availwritten;
%         1;
%     end
% end

% job done
status.jobdone = true;

end

