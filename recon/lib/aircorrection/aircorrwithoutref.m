function rawdata = aircorrwithoutref(rawdata, prmflow, KVmA, viewangle, aircorr)
% step 1 of air correction
% rawdata = aircorrwithoutref(rawdata, prmflow, KVmA, aircorr);

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
Nview_red = size(rawdata, 2);

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
    airKVmA = [airKVmA airKVmA(:, 1)];
else
    airKVmA = zeros(Nfocal, Nsect+1);
end

% interpolation index and weight
retangle = mod(viewangle - aircorr.firstangle, pi*2);
intp_index = floor(retangle./sectangle);
intp_alpha = retangle./sectangle - intp_index;
intp_index = intp_index + 1;

% KVmA
logKVmA = -log2(KVmA);

% corr rawdata with air
for ifocal = 1:Nfocal
    % indexes
    viewindex = ifocal:Nfocal:Nview_red;
    % un-matching focal spots' number warning
    ifocal_mod = mod((ifocal-1), double(aircorr.focalnumber)) + 1;
    airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal_mod-1);
    % rawdata
    rawdata(:, viewindex) = ...
        rawdata(:, viewindex) - airmain(airindex, intp_index(viewindex)).*(1-intp_alpha(viewindex));
    rawdata(:, viewindex) = ...
        rawdata(:, viewindex) - airmain(airindex, intp_index(viewindex)+1).*intp_alpha(viewindex);
    % KVmA
    logKVmA(viewindex) = ...
        logKVmA(viewindex) - airKVmA(ifocal_mod, intp_index(viewindex)).*(1-intp_alpha(viewindex)) ...
        - airKVmA(ifocal_mod, intp_index(viewindex)+1).*intp_alpha(viewindex);
end
1;
rawdata = rawdata - logKVmA;

% % default refblock
% Nref = prmflow.raw.air.refnumber;
% dataflow.rawhead.refblock = false(Nref, size(rawdata, 2));

end