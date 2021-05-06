function [dataflow, prmflow, status] = reconnode_Radialrebin(dataflow, prmflow, status)
% recon node, Axial rebin 
% [dataflow, prmflow, status] = reconnode_Radialrebin(dataflow, prmflow, status);

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
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nshot = prmflow.recon.Nshot;
% I know Nview=Nviewprot*Nshot, for axial

% rebin is prepared
rebin = prmflow.rebin;

% reshape rawdata 
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

% QDO reorder
if rebin.isQDO
    % reshape rawdata 
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice*Nview);
    % tmp
    A_QDO = zeros(rebin.Npixel, Nslice*Nview/2, 'single');
    % I know in QDO the rebin.Npixel~=Npixel
    Savail = ~isnan(rebin.QDOorder);
    % QDO phase1
    index_s1 = (1:Nslice*Nviewprot/2)' + (0:Nshot-1).*Nslice*Nviewprot;
    A_QDO(rebin.QDOorder(Savail(:,1), 1), :) = dataflow.rawdata(Savail(:,1), index_s1(:));
    % QDO phase2
    index_s2 = (Nslice*Nviewprot/2+1:Nslice*Nviewprot)' + (0:Nshot-1).*Nslice*Nviewprot;
    A_QDO(rebin.QDOorder(Savail(:,2), 2), :) = dataflow.rawdata(Savail(:,2), index_s2(:));
    % replace the rawdata
    dataflow.rawdata = reshape(A_QDO, rebin.Npixel*Nslice, Nview/2);
end

% radial rebin
dataflow.rawdata = dataflow.rawdata(rebin.radialindex(:), :).*(1-rebin.interalpha_rad(:)) + ...
    dataflow.rawdata(rebin.radialindex(:)+1,:).*rebin.interalpha_rad(:);

% prm to return
prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = rebin.Nviewprot;
prmflow.recon.Nview = rebin.Nviewprot*Nshot;
prmflow.recon.midchannel = rebin.midchannel;
prmflow.recon.delta_d = rebin.delta_d;
prmflow.recon.delta_z = rebin.delta_z;
prmflow.recon.SID = rebin.SID;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end