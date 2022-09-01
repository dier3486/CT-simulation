function [dataflow, prmflow, status] = reconnode_SlopeAzirebin(dataflow, prmflow, status)
% recon node, Axial slope-Azi rebin 
% [dataflow, prmflow, status] = reconnode_SlopeAzirebin(dataflow, prmflow, status);
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

Nshot =  prmflow.recon.Nshot;
Nslice = prmflow.recon.Nslice;
rebin = prmflow.rebin;
Nviewprot = single(rebin.Nviewprot);
Nreb = rebin.Nreb;
% I know it was Nviewprot/Nfocal.
idealphi = rebin.idealphi;
delta_view = single(pi*2/Nviewprot);
isGPU = ~isempty(status.GPUinfo);
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, []);
% prepare interp
if isGPU
    fAzi0 = gpuArray(-idealphi(:)./delta_view);
    fAzi = zeros(Nreb, Nviewprot, 'single', 'gpuArray');
    viewindex = gpuArray(single(1:Nviewprot+1));
    pixelindex = gpuArray(single(1:Nreb)');
else
    fAzi0 = -idealphi(:)./delta_view;
    viewindex = single(1:Nviewprot+1);
    pixelindex = single(1:Nreb)';
end

% viewangle
viewangle = reshape(dataflow.rawhead.viewangle, Nviewprot, Nshot);
startviewangle = viewangle(1,:);

% startvindex = mod(gather(max(floor(-fAzi))+1), Nviewprot) + 1;
% fAzi = mod(fAzi+viewindex([startvindex:Nviewprot 1:startvindex-1])-1, Nviewprot) + 1;

for ishot = 1:Nshot
    vindex = (1:Nviewprot)+(ishot-1)*Nviewprot;
    fAzi = fAzi0 + viewindex(1:end-1) + ((startviewangle(1)-startviewangle(ishot))/delta_view - 1);
    fAzi = mod(fAzi-1, Nviewprot) + 1;
    for islice = 1:Nslice
        % get data per slice per shot
        if isGPU
            data_islice = zeros(Nreb, Nviewprot+1, class(dataflow.rawdata), 'gpuArray');
        else
            data_islice = zeros(Nreb, Nviewprot+1, class(dataflow.rawdata));
        end
        data_islice(:, 1:Nviewprot) = squeeze(dataflow.rawdata(:, islice, vindex));
        % boundary
        data_islice(:, end) = data_islice(:, 1);
        % interp
        dataflow.rawdata(:, islice, vindex) = reshape(gather(interp2(viewindex, pixelindex, data_islice, ...
            fAzi, repmat(pixelindex, 1, Nviewprot), 'linear')), Nreb, 1, Nviewprot);
    end
    
end

% replace viewangle (?)
dataflow.rawhead.viewangle = repmat((0:Nviewprot-1).*(pi*2/Nviewprot)+startviewangle(1), 1, Nshot);

% recon coefficients
prmflow.recon.startviewangle = repmat(startviewangle(1), 1, Nshot);
prmflow.recon.delta_view = delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end