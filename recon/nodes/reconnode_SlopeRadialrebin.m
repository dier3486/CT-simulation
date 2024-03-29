function [dataflow, prmflow, status] = reconnode_SlopeRadialrebin(dataflow, prmflow, status)
% recon node, Axial slope-radial rebin 
% [dataflow, prmflow, status] = reconnode_SlopeRadialrebin(dataflow, prmflow, status);
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

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [prmflow, status] = reconnode_sloperebinprepare(prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

rebin = prmflow.rebin;
Nview =  rebin.Nview;
Nslice = rebin.Nslice;
Npixel = rebin.Npixel;
Nfocal = rebin.Nfocal;
% Nviewprot = prmflow.recon.Nviewprot;
viewblock = rebin.viewblock;
viewblock_focal = viewblock/Nfocal;

Nreb = rebin.Nreb;
isGPU = ~isempty(status.GPUinfo);

% get view angle (DFS)
% dataflow.rawhead.viewangle = dataflow.rawhead.viewangle(1:Nfocal:end) + (pi*2/Nviewprot)*(Nfocal-1)/2;
viewangle = dataflow.rawhead.viewangle(1:Nfocal:end) + rebin.DFSviewshift;

dataout = zeros(Nreb, Nslice, Nview/Nfocal, class(dataflow.rawdata));
if isGPU
    % GPU
    idealphi = gpuArray(rebin.idealphi);
    pixelindex = gpuArray(single(1:Npixel*Nfocal));
    faninterpkern = gpuArray(rebin.faninterpkern);
    viewangle = gpuArray(viewangle);
    Yshift = gpuArray(rebin.Yshift);
    dfanun1 = gpuArray(1/rebin.dfan);
    midU = gpuArray(rebin.midchannel);
    sliceindex = gpuArray(single(1:Nslice));
    viewindex = repmat(reshape(gpuArray(single(1:viewblock_focal)), 1, 1, []), Nreb, Nslice, 1);
    Zgrid = repmat(gpuArray(rebin.Zgrid), 1, 1,viewblock_focal);
    Xgrid = repmat(gpuArray(single(1:Nreb)'), 1, Nslice, viewblock_focal);
else
    idealphi = rebin.idealphi;
    pixelindex = single(1:Npixel*Nfocal);
    faninterpkern = rebin.faninterpkern;
    Yshift = rebin.Yshift;
    dfanun1 = 1/rebin.dfan;
    midU = rebin.midchannel;
    sliceindex = single(1:Nslice);
    viewindex = repmat(reshape(single(1:viewblock_focal), 1, 1, []), Nreb, Nslice, 1);
    Zgrid = repmat(rebin.Zgrid, 1, 1, viewblock_focal);
    Xgrid = repmat(single(1:Nreb)', 1, Nslice, viewblock_focal);
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
    % interp to equal fanangle
    for islice = 1:Nslice
        data_blk(:, islice, :) = reshape(interp1(pixelindex, reshape(data_blk(:, islice, :), Npixel*Nfocal, viewblock_focal),...
            faninterpkern(:, islice), 'linear', 'extrap'), Npixel*Nfocal, 1, viewblock_focal);
    end
    % interp to ideal fanangle (equal distance)
    phi3 = idealphi + reshape(asin(reshape(sin(viewangle(vindex_out)+idealphi), [], 1)*Yshift), Nreb, []);
    phi3 = phi3.*dfanun1 + midU;
    phi3 = permute(reshape(phi3, Nreb, viewblock_focal, Nslice), [1 3 2]);
    data2_blk = interp3(data_blk, repmat(sliceindex, Nreb, 1, viewblock_focal), phi3, viewindex, 'linear', 0);
    % interp Z
    dataout(:,:,vindex_out) = gather(interp3(data2_blk, Zgrid, Xgrid, viewindex, 'linear', 0));
end

% return data
dataflow.rawdata = dataout;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end