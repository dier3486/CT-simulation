function [dataflow, prmflow, status] = reconnode_FanRadialrebin(dataflow, prmflow, status)
% recon node, fan-radial rebin, the radial rebin before Azi
% [dataflow, prmflow, status] = reconnode_FanRadialrebin(dataflow, prmflow, status);
% Support (X)DFS

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

Nview =  prmflow.recon.Nview;
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nviewprot = prmflow.recon.Nviewprot;
% I know this Nviewprot is before rebin, which is not equal with rebin.Nviewprot
Nfocal = prmflow.recon.Nfocal;
rebin = prmflow.rebin;
viewblock = rebin.viewblock;
viewblock_focal = viewblock/Nfocal;

Nreb = rebin.Nreb;
isGPU = ~isempty(status.GPUinfo);

% view angle (DFS)
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(1:Nfocal:end) + (pi*2/Nviewprot)*(Nfocal-1)/2;

dataout = zeros(Nreb, Nslice, Nview/Nfocal, class(dataflow.rawdata));
if isGPU
    % GPU
    pixelindex = gpuArray(single(1:Npixel*Nfocal));
    faninterpkern = gpuArray(rebin.faninterpkern);
else
    pixelindex = single(1:Npixel*Nfocal);
    faninterpkern = rebin.faninterpkern;
end

dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

% DFS
if Nfocal>1
    for ifocal = 1:Nfocal
        interpshift = floor(rebin.DFSviewinterp(ifocal));
        interpalpha = rebin.DFSviewinterp(ifocal) - interpshift;
        
        index0 = ifocal:Nfocal:Nview;
        index1 = mod((1:Nview/Nfocal)+interpshift-1, Nview/Nfocal)+1;
        index1 = index0(index1);
        index2 = mod((1:Nview/Nfocal)+interpshift, Nview/Nfocal)+1;
        index2 = index0(index2);
        dataflow.rawdata(:, index0) = dataflow.rawdata(:, index1).*(1-interpalpha) + dataflow.rawdata(:, index2).*interpalpha;
    end
end

for iblk = 1:Nview/viewblock
    vindex_iblk = (1:viewblock)+(iblk-1)*viewblock;
    vindex_out = (1:viewblock_focal)+(iblk-1)*viewblock_focal;
    if isGPU
        data_blk = gpuArray(dataflow.rawdata(:, vindex_iblk));
    else
        data_blk = dataflow.rawdata(:, vindex_iblk);
    end
    % reorder, focal first
    data_blk = reshape(permute(reshape(data_blk, Npixel*Nslice, Nfocal, []), [2, 1, 3]), Npixel*Nfocal, Nslice, []);
    % interp 
    for islice = 1:Nslice
        data_tmp = reshape(interp1(pixelindex, reshape(data_blk(:, islice, :), Npixel*Nfocal, viewblock_focal),...
            faninterpkern(:, islice), 'linear', 'extrap'), Nreb*Nfocal, 1, viewblock_focal);
        dataout(:,islice,vindex_out) = gather(data_tmp);
    end
end

% return data
dataflow.rawdata = dataout;

% prm
prmflow.recon.Npixel = rebin.Nreb;
prmflow.recon.Nviewprot = rebin.Nviewprot;
prmflow.recon.Nview = Nview/Nfocal;
prmflow.recon.midchannel = rebin.midU_phi;
prmflow.recon.delta_d = rebin.delta_d;
prmflow.recon.delta_z = rebin.delta_z;
prmflow.recon.SID = rebin.SID;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end