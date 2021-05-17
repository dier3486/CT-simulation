function [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status)
% recon node, Azi rebin for axial
% [dataflow, prmflow, status] = reconnode_Azirebin(dataflow, prmflow, status);
% fly-focal is not supported yet

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

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nshot = prmflow.recon.Nshot;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nfocal = prmflow.recon.Nfocal;
% I know Nview=Nviewprot*Nshot, for axial
% to support DFS
Npixel_focal = Npixel*Nfocal;
Nviewprot_focal = Nviewprot/Nfocal;
Nview_focal = Nview/Nfocal;
% I know Nview_focal = Nview/Nfocal = Nviewprot_focal*Nshot, which is new view number after DFS Azirebin

% rebin is prepared
rebin = prmflow.rebin;

% I know for DFS the Npixel_rebin=Npixel*2, Nviewprot_rebin=Nviewprot*2.

% Azi rebin
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% debug
dataflow.rawdata(:,2:2:end) = 0;
for ishot = 1:Nshot
    % ini A (rawdata to be replace)
    A = zeros(Npixel_focal*Nslice, Nviewprot_focal, 'single');
    % views' index of input and output
    viewindex = (1:Nviewprot) + (ishot-1)*Nviewprot;
    viewindex_focal = (1:Nviewprot_focal) + (ishot-1)*Nviewprot_focal;
    % do interp
    A(rebin.vindex1_azi) = reshape(dataflow.rawdata(:, viewindex), [], Nviewprot_focal).*(1-rebin.interalpha_azi);
    A(rebin.vindex2_azi) = A(rebin.vindex2_azi) + ...
                           reshape(dataflow.rawdata(:, viewindex), [], Nviewprot_focal).*rebin.interalpha_azi;
    % start angle for first rebin view, replace the rawdata
    dataflow.rawdata(:, viewindex) = ...
        reshape([A(:, (rebin.startvindex : Nviewprot_focal)) A(:, (1 : rebin.startvindex-1))], Npixel*Nslice, Nviewprot);
    % I know the A is in shape (Npixel_rebin*Nslice, Nviewprot_rebin) but reshpaed to (Npixel*Nviewprot, Nviewprot) to fit for
    % dataflow.rawdata, which will be reshaped back after looping the shots
    
    % ref blocked
    % to find out the views whose references are blocked after rebin
    if isfield(dataflow.rawhead, 'refblock') && any(dataflow.rawhead.refblock(:))
        % B is the refblock to be replace
        B = false(2, Nviewprot);
        Bindex_shift = (0:Nfocal-1)'.*Nviewprot_focal;
        % ref1
        index_ref1 = (0:Nfocal-1).*(Npixel*Nslice) + 1;
        v1_ref1 = ceil(rebin.vindex1_azi(index_ref1, :) ./ (Npixel*Nslice*Nfocal));
        v2_ref1 = ceil(rebin.vindex2_azi(index_ref1, :) ./ (Npixel*Nslice*Nfocal));
        B(1, v1_ref1 + Bindex_shift) = dataflow.rawhead.refblock(1, viewindex);
        B(1, v2_ref1 + Bindex_shift) = B(1, v2_ref1 + Bindex_shift) | dataflow.rawhead.refblock(1, viewindex);
        % ref2
        index_ref2 = (1:Nfocal).*(Npixel*Nslice);
        v1_ref2 = ceil(rebin.vindex1_azi(index_ref2, :) ./ (Npixel*Nslice*Nfocal));
        v2_ref2 = ceil(rebin.vindex2_azi(index_ref2, :) ./ (Npixel*Nslice*Nfocal));
        B(2, v1_ref2 + Bindex_shift) = dataflow.rawhead.refblock(2, viewindex);
        B(2, v2_ref2 + Bindex_shift) = B(2, v2_ref2 + Bindex_shift) | dataflow.rawhead.refblock(2, viewindex);
        % do any for DFS
        B = any(reshape(B, 2, Nviewprot_focal, Nfocal), 3);
        % start angle
        dataflow.rawhead.refblock(:, viewindex_focal) = ...
            [B(:, (rebin.startvindex : Nviewprot_focal)) B(:, (1 : rebin.startvindex-1))];
    end
end
% reshape A for DFS
dataflow.rawdata = reshape(dataflow.rawdata, Npixel_focal*Nslice, Nview_focal);
% cut B for DFS
dataflow.rawhead.refblock = dataflow.rawhead.refblock(:, 1:Nview_focal);

% viewangle
viewangle = reshape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
viewangle = viewangle(1:Nfocal:end, :);
startviewangle = viewangle(rebin.startvindex, :);
dataflow.rawhead.viewangle = [viewangle(rebin.startvindex : Nviewprot_focal, :); viewangle(1 : rebin.startvindex-1, :)];
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(:)';

% prm to return
prmflow.recon.Nview = Nview_focal;
prmflow.recon.Npixel = Npixel_focal;
prmflow.recon.Nviewprot = Nviewprot_focal;
prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = rebin.delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end