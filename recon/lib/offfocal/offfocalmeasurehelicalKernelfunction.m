function offraw = offfocalmeasurehelicalKernelfunction(dataflow, prmflow, status, buffer)
% offfocal correction step2

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

% range in reading the dataflow.rawdata
if pipeline_onoff
    % new offraw write view index (sparsed) 
    offWriteViewindex = buffer.offReadViewindex + buffer.offWritePoint - buffer.offReadPoint;
    % input data start view index
    datastartView = (offWriteViewindex - extraview(2) - 1)*viewsparse + 1;
    % input data start point (whose view index is datastartView)
    datastartPoint = buffer.AvailPoint - buffer.AvailViewindex + datastartView;
    % move the start point if less than AvailPoint
    datastartPoint = max(datastartPoint, buffer.AvailPoint + 1);    % never to read the data before AvailPoint
    % the datastartView will be moved to
    datastartView = buffer.AvailViewindex - buffer.AvailPoint + datastartPoint;
    % input data end point
    dataendPoint = buffer.WritePoint - 1;
    dataendView = buffer.AvailViewindex - buffer.AvailPoint + dataendPoint;
    % new raw1 end view index (sparsed)
    if dataendView < Nview
        offEndViewindex = floor((dataendView + extraview(1) - 1) / viewsparse) + 1;
    else
        offEndViewindex = prmoff.offendview;
    end
else
    datastartPoint = 1;
    dataendPoint = Nview;
end
% view number to be used
Nrenew = dataendPoint - datastartPoint + 1;

% I know this was done,
% dataflow.rawdata = 2.^(-dataflow.rawdata);

% merge slices
if slicemerge>1
    if ~slicezebra
        Aoff = squeeze(mean(reshape(dataflow.rawdata(:, datastartPoint:dataendPoint), ...
            Npixel, slicemerge, Nslicemerge, Nrenew), 2));
    else
        Aoff = reshape(mean(reshape(dataflow.rawdata(:, datastartPoint:dataendPoint), ...
            Npixel, 2, slicemerge, Nslicemerge/2, Nrenew), 3), Npixel, Nslicemerge, Nrenew);
    end
else
    Aoff = reshape(dataflow.rawdata(:, datastartPoint:dataendPoint), Npixel, Nslice, Nrenew);
end

% range in writing the 'offraw' in off-focal space
if pipeline_onoff
    % new offraw start reletive view index
    offstartView = offWriteViewindex;
    offendView = offEndViewindex;
    index_voff = ((offstartView : offendView) - 1).*viewsparse + 1 - datastartView + 1;
    Nviewoff = max(offendView - offstartView + 1, 0);
else
    offstartView = prmoff.offstartview;
    offendView = prmoff.offendview;
    index_voff = ((offstartView : offendView) - 1).*viewsparse + 1;
    Nviewoff = prmoff.Nviewoff;
end

% interp Aoff to off-focal measure space
offraw = struct();
offraw.rawdata = zeros(Noffsample, Nslicemerge, Nviewoff, 'single');
for islice = 1:Nslicemerge
    Df = index_voff - prmoff.rawinterp2phi;
    Df(Df<1) = 1;   Df(Df>Nrenew) = Nrenew;
    offraw.rawdata(:, islice, :) = interp2(squeeze(Aoff(:, islice, :)) , Df, repmat(prmoff.rawinterp2t(:, islice), 1, Nviewoff), 'linear', 0);
end

% conv
offraw.rawdata = ifft(fft(reshape(offraw.rawdata, Noffsample, []), Noffsample).*prmoff.offkernel);
offraw.rawdata = reshape(offraw.rawdata, Noffsample*Nslicemerge, Nviewoff);

end