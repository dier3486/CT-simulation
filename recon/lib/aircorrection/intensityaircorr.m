function dataflow = intensityaircorr(dataflow, prmflow, aircorr)
% intensity air correction, used in calibration
% dataflow = intensityaircorr(dataflow, prmflow, aircorr);

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

% parameters to use in prmflow
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
% Nview = prmflow.raw.Nview;
Nfocal = prmflow.raw.Nfocal;
Nview_red = size(dataflow.rawdata, 2);

% parameters in corr
Nsect = single(aircorr.Nsection);

% angles of the air table
sectangle = (pi*2/Nsect);
% airangle = (-1:Nsect).*(pi*2/Nsect);
% airmain & airref
aircorr.main = reshape(aircorr.main, [], Nsect);
airmain = [aircorr.main aircorr.main(:,1)];
if isfield(aircorr, 'referenceKVmA')
    airKVmA = reshape(aircorr.referenceKVmA, [], Nsect);
    if size(airKVmA, 1) < Nfocal
        % This should be a bug in air calibration table
        airKVmA = repmat(airKVmA, Nfocal, 1);
    end
    airKVmA = [airKVmA airKVmA(:, 1)];
else
    airKVmA = nan(Nfocal, Nsect+1);
end

% interpolation index and weight
retangle = mod(dataflow.rawhead.viewangle - aircorr.firstangle, pi*2);
intp_index = floor(retangle./sectangle);
intp_alpha = retangle./sectangle - intp_index;
intp_index = intp_index + 1;

% KVmA
KVmA = dataflow.rawhead.KV.*dataflow.rawhead.mA;

% offset correction
dataflow.rawdata = dataflow.rawdata - mean(dataflow.offset.rawdata, 2);

% Integration Time
inttime = mean(single(dataflow.rawhead.Integration_Time));

% corr rawdata with air
for ifocal = 1:Nfocal
    % rawdata
    viewindex = ifocal:Nfocal:Nview_red;
    airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal-1);
    Iair = 2.^(-airmain(airindex, intp_index(viewindex)).*(1-intp_alpha(viewindex)) - ...
        airmain(airindex, intp_index(viewindex)+1).*intp_alpha(viewindex)).*inttime;
    dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex)./Iair;

    % KVmA
    KVmA(viewindex) = ...
        (airKVmA(ifocal, intp_index(viewindex)).*(1-intp_alpha(viewindex)) ...
        + airKVmA(ifocal, intp_index(viewindex)+1).*intp_alpha(viewindex))./KVmA(viewindex);
end
dataflow.rawdata = dataflow.rawdata./KVmA;

% default refblock
Nref = prmflow.raw.air.refnumber;
dataflow.rawhead.refblock = false(Nref, size(dataflow.rawdata, 2));

end