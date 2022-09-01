function [dataflow, prmflow, status] = reconnode_HelicalAzirebin(dataflow, prmflow, status)
% recon node, Helical Azi rebin 
% [dataflow, prmflow, status] = reconnode_HelicalAzirebin(dataflow, prmflow, status);
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

Nview = prmflow.recon.Nview;
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
    fAzi = fAzi0 + gpuArray(single(0:Nview-1));
    viewindex = gpuArray(single(1:Nview));
    pixelindex = gpuArray(single(1:Nreb)');
else
    fAzi0 = -idealphi(:)./delta_view;
    fAzi = fAzi0 + single(0:Nview-1);
    viewindex = single(1:Nview);
    pixelindex = single(1:Nreb)';
end

% start viewangle
startviewangle = dataflow.rawhead.viewangle(1);

% startvindex = mod(gather(max(floor(-fAzi))+1), Nviewprot) + 1;
% fAzi = mod(fAzi+viewindex([startvindex:Nviewprot 1:startvindex-1])-1, Nviewprot) + 1;

for islice = 1:Nslice
    % get data per slice per shot
    if isGPU
        data_islice = gpuArray(squeeze(dataflow.rawdata(:, islice, :)));
    else
        data_islice = squeeze(dataflow.rawdata(:, islice, :));
    end 
    % interp
%     dataflow.rawdata(:, islice, :) = reshape(fillmissing(gather(interp2(viewindex, pixelindex, data_islice, ...
%         fAzi, repmat(pixelindex, 1, Nview), 'linear')), 'nearest', 2), Nreb, 1, Nview);
    dataflow.rawdata(:, islice, :) = reshape(gather(interp2(viewindex, pixelindex, data_islice, ...
        fAzi, repmat(pixelindex, 1, Nview), 'linear', 0)), Nreb, 1, Nview);

end

% recon coefficients
prmflow.recon.startviewangle = startviewangle;
prmflow.recon.delta_view = delta_view;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end