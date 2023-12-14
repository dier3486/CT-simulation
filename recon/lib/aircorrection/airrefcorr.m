function [dataflow, reflast] = airrefcorr(dataflow, prmflow, altref, startrefviewindex, reflast, sliceindependent)
% step 2 of air correction
% [dataflow, reflast] = airrefcorr(dataflow, prmflow, altref, startrefviewindex, reflast, sliceindependent);

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

Npixel = prmflow.raw.Npixel;
% Nslice = prmflow.raw.Nslice;
% Nview = prmflow.raw.Nview;
Nfocal = prmflow.raw.Nfocal;
% Nshot = prmflow.raw.Nshot;
blockwindow = prmflow.raw.air.blockwindow;

shots = unique(dataflow.rawhead.Shot_Number);

for ishots = length(shots)
    % clear when ishots>1
    if ishots > 1
        startrefviewindex = 1;
        reflast = cell(1, Nfocal);
    end

    % view index of current shot
    viewindex = find(dataflow.rawhead.Shot_Number == shots(ishots));

    dataflow.rawhead.refblock(:, viewindex) = isrefblocked(dataflow.rawhead.refblock(:, viewindex), prmflow, ...
        altref.referr(:, viewindex), startrefviewindex);

    % view delay for axial
    if strcmpi(prmflow.raw.scan, 'axial') && isempty(startrefviewindex)
        viewindex = [viewindex(blockwindow+1 : end) viewindex(1 : blockwindow)];
    end

    % corr rawdata with ref
    for ifocal = 1:Nfocal
        viewindex_ifc = viewindex(ifocal:Nfocal:end);
        [rawref_ifc, reflast{ifocal}] = referencecorr(altref.rawref(:, viewindex_ifc), ...
            dataflow.rawhead.refblock(:, viewindex_ifc), ...
            altref.KVmA(:, viewindex_ifc), reflast{ifocal});
        if sliceindependent
            dataflow.rawdata(:, viewindex_ifc) = dataflow.rawdata(:, viewindex_ifc) - repelem(rawref_ifc, Npixel, 1);
        else
            dataflow.rawdata(:, viewindex_ifc) = dataflow.rawdata(:, viewindex_ifc) - rawref_ifc;
        end
    end

end

end