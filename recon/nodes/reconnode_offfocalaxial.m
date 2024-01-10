function [dataflow, prmflow, status] = reconnode_offfocalaxial(dataflow, prmflow, status)
% recon node, off-focal correction for axial
% [dataflow, prmflow, status] = reconnode_offfocalaxial(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_offfocalprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters to use in prmflow
Nshot = prmflow.raw.Nshot;
Nview = prmflow.raw.Nview;
Nslice = prmflow.raw.Nslice;
Npixel = prmflow.raw.Npixel;
Nviewprot = prmflow.raw.Nviewprot;
% prepared parameters
prmoff = prmflow.raw.offfocal;
slicemerge = prmoff.slicemerge;
Nslicemerge = Nslice/slicemerge;
Noffsample = prmoff.Noffsample;
Nviewoff = prmoff.Nviewoff;
slicezebra = prmoff.slicezebra;

% exp
dataflow.rawdata = 2.^(-dataflow.rawdata);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nviewprot, Nshot);

Aoff = zeros(Npixel, Nslicemerge, Nviewprot+1, 'single');
for ishot = 1:Nshot
    % slice merge
    if slicemerge>1
        if ~slicezebra
            Aoff(:, :, 1:end-1) = squeeze(mean( ...
                reshape(dataflow.rawdata(:, :, ishot), Npixel, slicemerge, Nslicemerge, Nviewprot), 2));
        else
            Aoff(:, :, 1:end-1) = reshape(mean( ...
                reshape(dataflow.rawdata(:, :, ishot), Npixel, 2, slicemerge, Nslicemerge/2, Nviewprot), 3), ...
                Npixel, Nslicemerge, Nviewprot);
        end
    else
        Aoff(:, :, 1:end-1) = reshape(dataflow.rawdata(:, :, ishot), Npixel, Nslice, Nviewprot);
    end
    Aoff(:, :, end) = Aoff(:, :, 1);

    % off-focal view index
    index_voff = (0 : prmoff.offendview-prmoff.offstartview).*prmoff.viewsparse + prmoff.offstartview;
    % raw view index
    index_vraw = single(0:Nviewprot-1)./prmoff.viewsparse + 1;

    % interp raw0 to off-focal measure space
    raw1 = zeros(Noffsample, Nslicemerge, Nviewoff+1, 'single');
    for islice = 1:Nslicemerge
        Df = mod(index_voff - prmoff.rawinterp2phi - 1, Nviewprot) + 1;
        raw1(:, islice, 1:Nviewoff) = interp2(squeeze(Aoff(:, islice, :)) , Df, ...
            repmat(prmoff.rawinterp2t(:, islice), 1, Nviewoff), 'linear', 0);
    end
    raw1(:, :, end) =  raw1(:, :, 1);
    1;
    % conv
    raw1 = ifft(fft(reshape(raw1, Noffsample, []), Noffsample).*prmoff.offkernel);
    raw1 = reshape(raw1, Noffsample, Nslicemerge, []);

    % interp raw1 back to raw space
    Afix = zeros(Npixel, Nslicemerge, Nviewprot, 'single');
    for islice = 1:Nslicemerge
        Dfb = mod(index_vraw + prmoff.tinterp2phi(:, islice) - 1, Nviewoff) + 1;
        Afix(:, islice, :) = interp2(squeeze(raw1(:, islice, :)) , Dfb, ...
            repmat(prmoff.tinterp2raw(:, islice), 1, Nviewprot), 'linear', 0);
    end
    1;
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
    Afix = reshape(permute(reshape(Afix, Nslice, Npixel, Nviewprot), [2 1 3]), Nslice*Npixel, Nviewprot);

    % add to rawdata
    dataflow.rawdata(:, :, ishot) = dataflow.rawdata(:, :, ishot) - Afix;
end

% minimum boundary
minintensity = prmflow.raw.offfocal.minintensity;
dataflow.rawdata(dataflow.rawdata<minintensity) = minintensity;
% log2
dataflow.rawdata = -log2(dataflow.rawdata);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end