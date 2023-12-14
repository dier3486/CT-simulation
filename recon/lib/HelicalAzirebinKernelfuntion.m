function [dataflow, prmflow, status] = HelicalAzirebinKernelfuntion(dataflow, prmflow, status)
% recon node, helical azi-rebin kernel function
% [dataflow, prmflow, status] = HelicalAzirebinKernelfuntion(dataflow, prmflow, status);

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

% nodename
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% parameters from rebin
rebin = prmflow.rebin;

% Nextraview
if pipeline_onoff && isfield(rebin, 'Nextraview')
    Nextraview = rebin.Nextraview;
else
    Nextraview = [0 0];
end

Nslice = rebin.Nslice;
% Nfocal = rebin.Nfocal;
% Nview = rebin.Nview / Nfocal;
% Nviewprot = rebin.Nviewprot / Nfocal;
Nreb = rebin.Nreb;
idealphi = rebin.idealphi;
delta_view = rebin.delta_view;
isGPU = ~isempty(status.GPUinfo);
Nview = size(dataflow.rawdata, 2);

viewindex_out = single(1-Nextraview(1) : Nview+Nextraview(2));
Nview_out = Nview + Nextraview(1) + Nextraview(2);

% viewangle shift pi/2, be used to
dataflow.rawhead.viewangle = dataflow.rawhead.viewangle + pi/2;
% I know the rawhead.viewangle has been replaced in reconnode_FanRadialrebin

% view extra
dataflow = dataviewextend(dataflow, Nextraview);    

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Nreb, Nslice, []);

% fAzi
fAzi0 = -idealphi(:)./delta_view;
% low accuracy interpolation
if isfield(rebin, 'lowaccuracyshift')
    fAzi0 = round(fAzi0.*2^rebin.lowaccuracyshift).*2^(-rebin.lowaccuracyshift);
end

% gpu array
if isGPU
    fAzi = gpuArray(fAzi0) + gpuArray(viewindex_out);
    viewindex = gpuArray(single(0:Nview+1));
    pixelindex = gpuArray(single(1:Nreb)');
    data_islice = zeros(Nreb, Nview+2, 'single', 'gpuArray');
else
    fAzi = fAzi0 + viewindex_out;
    viewindex = single(0:Nview+1);
    pixelindex = single(1:Nreb)';
    data_islice = zeros(Nreb, Nview+2, 'single');
end

for islice = 1:Nslice
    % get data per slice per shot
    if isGPU
        data_islice(:, 2:end-1) = gpuArray(squeeze(dataflow.rawdata(:, islice, Nextraview(1)+1 : Nextraview(1)+Nview)));
    else
        data_islice(:, 2:end-1) = squeeze(dataflow.rawdata(:, islice, Nextraview(1)+1 : Nextraview(1)+Nview));
    end 
    1;
    % interp
%     dataflow.rawdata(:, islice, :) = reshape(fillmissing(gather(interp2(viewindex, pixelindex, data_islice, ...
%         fAzi, repmat(pixelindex, 1, Nview), 'linear')), 'nearest', 2), Nreb, 1, Nview);
    dataflow.rawdata(:, islice, :) = reshape(gather(interp2(viewindex, pixelindex, data_islice, ...
        fAzi, repmat(pixelindex, 1, Nview_out), 'linear', 0)), Nreb, 1, Nview_out);

end
1;
% reshape
dataflow.rawdata = reshape(dataflow.rawdata, [], Nview_out);


end

function dataflow = dataviewextend(dataflow, Nextraview)

if all(Nextraview==0)
    return;
end

NextraMax = max(Nextraview, [0 0]);
Nextra0 = min(Nextraview, [0 0]);

size1 = size(dataflow.rawdata, 1);
dataflow.rawdata = cat(2, zeros(size1, NextraMax(1), 'single'), dataflow.rawdata(:, 1-Nextra0(1) : end+Nextra0(2)), ...
    zeros(size1, NextraMax(2), 'single'));

headfields = fieldnames(dataflow.rawhead);
for ii = 1:length(headfields)
    size_ii = size(dataflow.rawhead.(headfields{ii}), 1);
    dataflow.rawhead.(headfields{ii}) = cat(2, zeros(size_ii, NextraMax(1), 'like', dataflow.rawhead.(headfields{ii})), ...
        dataflow.rawhead.(headfields{ii})(:, 1-Nextra0(1) : end+Nextra0(2)), ...
        zeros(size_ii, NextraMax(2), 'like', dataflow.rawhead.(headfields{ii})));
end

end
