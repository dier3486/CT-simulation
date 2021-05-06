function [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status)
% recon node, crosstalk correction
% [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status);

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
Nview = prmflow.recon.Nview;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;

% parameters set in pipe
crossprm = prmflow.pipe.(status.nodename);
if isfield(crossprm, 'weight')
    weight = crossprm.weight;
else
    weight = 1.0;
end
if isfield(crossprm, 'istointensity')
    istointensity = crossprm.istointensity;
else
    % defualt, do crosstalk in intensity
    istointensity = true;
end

% calibration table
crscorr = prmflow.corrtable.(status.nodename);
crsorder = crscorr.order;
% DFS
if isfield(crscorr, 'focalnumber') && crscorr.focalnumber
    Nfocal = crscorr.focalnumber;
else
    Nfocal = 1;
end
crsval_org = reshape(crscorr.main, [], crsorder, Nfocal);

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% to intensity
if istointensity
    dataflow.rawdata = 2.^(-dataflow.rawdata);
end
% correct
for ifocal = 1:Nfocal
    viewindex = ifocal:Nfocal:Nview;
    crsval = crsval_org(:, :, ifocal);
    switch crsorder
        case 1
            % odd-symmetric style
            % the correction operator is a tridiagonal matrix [-crsval;  1; crsval];
            % rawfix
            rawfix = zeros(Npixel*Nslice, Nview/Nfocal);
            rawfix(2:end-1, :) = (dataflow.rawdata(3:end, viewindex) - dataflow.rawdata(1:end-2, viewindex)).*crsval(2:end-1);
            % add to rawdata
            dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) + rawfix.*weight;
        case 3
            % 0
            rawfix = dataflow.rawdata(:, viewindex).*crsval(:,2);
            % -1
            rawfix(2:end, :) = rawfix(2:end, :) + dataflow.rawdata(1:end-1, viewindex).*crsval(2:end, 1);
            % 1
            rawfix(1:end-1, :) = rawfix(1:end-1, :) + dataflow.rawdata(2:end, viewindex).*crsval(1:end-1, 3);
            % add to rawdata
            dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) + rawfix.*weight;
        case 5
            % 0
            rawfix = dataflow.rawdata(:, viewindex).*crsval(:,3);
            % -2
            rawfix(3:end, :) = rawfix(3:end, :) + dataflow.rawdata(1:end-2, viewindex).*crsval(3:end, 1);
            % -1
            rawfix(2:end, :) = rawfix(2:end, :) + dataflow.rawdata(1:end-1, viewindex).*crsval(2:end, 2);
            % 1
            rawfix(1:end-1, :) = rawfix(1:end-1, :) + dataflow.rawdata(2:end, viewindex).*crsval(1:end-1, 4);
            % 2
            rawfix(1:end-2, :) = rawfix(1:end-2, :) + dataflow.rawdata(3:end, viewindex).*crsval(1:end-2, 5);
            % add to rawdata
            dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) + rawfix.*weight;

        otherwise
            error('Crosstalk correction not support order %d.', crsorder);
    end
end

if istointensity
    % min cut
    minval = 2^-32;
    dataflow.rawdata(dataflow.rawdata<minval) = minval;
    % log2
    dataflow.rawdata = -log2(dataflow.rawdata);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end