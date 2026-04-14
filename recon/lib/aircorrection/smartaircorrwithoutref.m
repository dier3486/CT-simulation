function rawdata = smartaircorrwithoutref(rawdata, prmraw, viewangle, smartaircorr)
% step 1 of smartair correction
% rawdata = smartaircorrwithoutref(rawdata, prmflow, viewangle, smartaircorr);

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
Npixel = prmraw.Npixel;
Nslice = prmraw.Nslice;
% Nview = prmflow.raw.Nview;
Nfocal = prmraw.Nfocal;
Nview_red = size(rawdata, 2);

% parameters in corr
Nsect = single(smartaircorr.Nsection);

% angles of the air table
sectangle = (pi*2/Nsect);
% airangle = (-1:Nsect).*(pi*2/Nsect);
% airmain & airref
smartaircorr.main = reshape(smartaircorr.main, [], Nsect);
smartairmain = cast([smartaircorr.main smartaircorr.main(:,1)], 'like', rawdata);

% interpolation index and weight
retangle = mod(viewangle, pi*2);
intp_index = floor(retangle./sectangle);
intp_alpha = retangle./sectangle - intp_index;
intp_index = intp_index + 1;

% corr rawdata with air
for ifocal = 1:Nfocal
    % indexes
    viewindex = ifocal:Nfocal:Nview_red;
    % un-matching focal spots' number warning
    ifocal_mod = mod((ifocal-1), double(smartaircorr.focalnumber)) + 1;  % int
    airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal_mod-1);
    % rawdata
    rawdata(:, viewindex) = ...
        rawdata(:, viewindex) - smartairmain(airindex, intp_index(viewindex)).*(1-intp_alpha(viewindex));
    rawdata(:, viewindex) = ...
        rawdata(:, viewindex) - smartairmain(airindex, intp_index(viewindex)+1).*intp_alpha(viewindex);
end

end