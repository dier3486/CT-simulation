function [dataflow, Nrenew, offNremove] = offfocalfixhelicalKernelfunction(dataflow, prmflow, status, offraw, buffer)
% offfocal correction step3

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

% parameters to use in prmflow
Nview = prmflow.raw.Nview;
Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
% prepared parameters
prmoff = prmflow.raw.offfocal;
slicemerge = prmoff.slicemerge;
Nslicemerge = Nslice/slicemerge;
Noffsample = prmoff.Noffsample;

slicezebra = prmoff.slicezebra;
extraview = prmoff.extraview;
viewsparse = prmoff.viewsparse;

if pipeline_onoff
    % offraw end view index (sparsed) 
    offEndViewindex = buffer.offReadViewindex + buffer.offWritePoint - buffer.offReadPoint - 1;
    % output data end view index
    if offEndViewindex < prmoff.offendview
        dataendView = (offEndViewindex - extraview(2) - 1) * viewsparse + 1;
    else
        dataendView = Nview;
    end
    % offraw start point
    offstartPoint = buffer.offReadPoint;
    % offraw end point
    offendPoint = buffer.offWritePoint - 1;
    Nviewoff = offendPoint - offstartPoint + 1;
    % output data start point
    datastartPoint = buffer.AvailPoint + 1;
    % output data end point
    dataendPoint = buffer.AvailPoint - buffer.AvailViewindex + dataendView;
    % view index to interp to raw space
    index_vraw = (buffer.AvailViewindex : dataendView-1)./viewsparse - buffer.offReadViewindex + 2;
else
    Nviewoff = size(offraw.rawdata, 2);
    offstartPoint = 1;
    offendPoint = Nviewoff;
    datastartPoint = 1;
    dataendPoint = Nview;
    % view index to interp to raw space
    index_vraw = (0:Nview-1)./viewsparse - prmoff.offstartview + 2;
end
% the view number to be corrected, will be returned to move the related points
Nrenew = dataendPoint - datastartPoint + 1;

% I know the new buffer.AvailPoint shall be buffer.AvailPoint+Nrenew
% the new buffer.offReadPoint will be buffer.offReadPoint+offNrenew, by
if pipeline_onoff
    offNremove = floor(dataendView/viewsparse) + 1 + extraview(1) - buffer.offReadViewindex;
else
    offNremove = 0;
end

% reshape
offraw.rawdata = reshape(offraw.rawdata, Noffsample, Nslicemerge, []);

% interp raw1 back to raw space
Afix = zeros(Npixel, Nslicemerge, Nrenew, 'single');
for islice = 1:Nslicemerge
    Dfb = index_vraw + prmoff.tinterp2phi(:, islice);
    Dfb(Dfb<1) = 1;   Dfb(Dfb>Nviewoff) = Nviewoff;
    Afix(:, islice, :) = interp2(squeeze(offraw.rawdata(:, islice, offstartPoint:offendPoint)) , Dfb, ...
        repmat(prmoff.tinterp2raw(:, islice), 1, Nrenew), 'linear', 0);
end

% measure scale
Afix = Afix.*prmoff.Dphiscale;
Afix = real(Afix) + imag(Afix).*prmoff.Dphiscale_odd;

% permute Afix to move the slice to dim 1
Afix = reshape(permute(Afix, [2 1 3]), Nslicemerge, []);
% inverse slice merge
if slicemerge>1
    Afix = repelem(Afix, slicemerge, 1);
    if slicezebra
        Afix = reshape(permute(reshape(Afix, slicemerge, 2, []), [2 1 3]), Nslice, []);
    end
end

% Z cross
if ~slicezebra
    Afix = prmoff.crsMatrix * Afix;
else
    Afix(1:2:end, :) = prmoff.crsMatrix * Afix(1:2:end, :);
    Afix(2:2:end, :) = prmoff.crsMatrix * Afix(2:2:end, :);
end

% reshape
Afix = reshape(permute(reshape(Afix, Nslice, Npixel, Nrenew), [2 1 3]), Nslice*Npixel, Nrenew);

% fix to rawdata
Afix = dataflow.rawdata(:, datastartPoint: dataendPoint) - Afix;

% minimum boundary
minintensity = prmflow.raw.offfocal.minintensity;
Afix(Afix<minintensity) = minintensity;
% log2
dataflow.rawdata(:, datastartPoint: dataendPoint) = -log2(Afix);

end