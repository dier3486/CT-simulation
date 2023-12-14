function [dataflow, prmflow, status] = FanRadialrebinKernelfunction(dataflow, prmflow, status)
% recon node, fan-radial rebin kernel function
% [dataflow, prmflow, status] = FanRadialrebinKernelfunction(dataflow, prmflow, status);

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

rebin = prmflow.rebin;
% Nview =  rebin.Nview;
Nslice = rebin.Nslice;
Npixel = rebin.Npixel;
Nfocal = rebin.Nfocal;
% Nviewprot = rebin.Nviewprot;
% I know this Nviewprot is before rebin, which shall not equal to recon.Nviewprot
viewblock = rebin.viewblock;
viewblock_focal = viewblock/Nfocal;

Nview = size(dataflow.rawhead.viewangle, 2);
Nblk = ceil(Nview/viewblock);

Nreb = rebin.Nreb;
isGPU = ~isempty(status.GPUinfo);

% Nview = length(dataflow.rawhead.viewangle);
% view angle (DFS)
dataflow.rawhead = rawheadfanradialrebin(dataflow.rawhead, rebin);

dataout = zeros(Nreb, Nslice, Nview/Nfocal, class(dataflow.rawdata));
if isGPU
    % GPU
    pixelindex = gpuArray(single(1:Npixel*Nfocal));
    faninterpkern = gpuArray(rebin.faninterpkern);
else
    pixelindex = single(1:Npixel*Nfocal);
    faninterpkern = rebin.faninterpkern;
end

% dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

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

for iblk = 1: Nblk
    if iblk < Nblk
        vindex_iblk = (1:viewblock)+(iblk-1)*viewblock;
        vindex_out = (1:viewblock_focal)+(iblk-1)*viewblock_focal;
        viewblock_focal_iblk = viewblock_focal;
    else
        vindex_iblk = viewblock*(Nblk-1)+1 : Nview;
        vindex_out = viewblock_focal*(Nblk-1)+1 : Nview/Nfocal;
        viewblock_focal_iblk = Nview/Nfocal - viewblock_focal*(Nblk-1);
    end
    if isGPU
        data_blk = gpuArray(dataflow.rawdata(:, vindex_iblk));
    else
        data_blk = dataflow.rawdata(:, vindex_iblk);
    end
    % reorder, focal first
    data_blk = reshape(permute(reshape(data_blk, Npixel*Nslice, Nfocal, []), [2, 1, 3]), Npixel*Nfocal, Nslice, []);
    % interp 
    for islice = 1:Nslice
        data_tmp = reshape(interp1(pixelindex, reshape(data_blk(:, islice, :), Npixel*Nfocal, viewblock_focal_iblk),...
            faninterpkern(:, islice), 'linear', 'extrap'), Nreb, 1, viewblock_focal_iblk);
        dataout(:,islice,vindex_out) = gather(data_tmp);
    end
end

% return data
dataflow.rawdata = reshape(dataout, [], Nview/Nfocal);

end