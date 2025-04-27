function [dataflow, prmflow, status] = reconnode_idealwater(dataflow, prmflow, status)
% cali node, get projection of ideal water phantom by a real scan
% [dataflow, prmflow, status] = reconnode_idealwater(dataflow, prmflow,
% status);
% use to find the ideal target for BH and nonlinear calibration.

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

% parameters in pipe
caliprm = prmflow.pipe.(status.nodename);
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    idealwater = prmflow.corrtable.(status.nodename);
    idealwater = everything2single(idealwater, 'any', 'single');
    tooutputtcorr = false;
else
    tooutputtcorr = true;
end

if tooutputtcorr
    % get ideal water
%     [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
    [dataflow, prmflow, status] = reconnode_Rebin(dataflow, prmflow, status);

    [dataflow, prmflow, status] = reconnode_watergoback(dataflow, prmflow, status);
    [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status);
    % to corr
    dataflow.idealwatercorr = caliprmforcorr(prmflow, corrversion);
    dataflow.idealwatercorr.Nslice = prmflow.rebin.Nslice;
    dataflow.idealwatercorr.viewnumber = prmflow.rebin.Nview;
    dataflow.idealwatercorr.mainsize = size(dataflow.rawdata(:), 1);
    dataflow.idealwatercorr.viewangle = dataflow.rawhead.viewangle;
    dataflow.idealwatercorr.indexrange = dataflow.rawhead.index_range;
    dataflow.idealwatercorr.main = dataflow.rawdata;
else
    % get ideal water from corr
    % to intercept the corr and reshape
    Nslice = prmflow.raw.Nslice;
    sliceindex = (1:Nslice) + (idealwater.Nslice-Nslice)*2;
    idealwater.main = reshape(idealwater.main, idealwater.Npixel, idealwater.Nslice, idealwater.viewnumber);
    idealwater.main = reshape(idealwater.main(:, sliceindex, :), idealwater.Npixel*Nslice, idealwater.viewnumber);
    % view interp
    [~, indexmax] = max(idealwater.viewangle);
    indexmin = mod(indexmax, idealwater.viewnumber) + 1;
    index_sort = [indexmin:idealwater.viewnumber  1:indexmin-1];
    view0 = [idealwater.viewangle(indexmax)-pi*2; idealwater.viewangle(index_sort); idealwater.viewangle(indexmin)+pi*2];
    % interp prepare
    [interpindex, interpalpha] = interpprepare(view0, dataflow.rawhead.viewangle(:));
    % reverse index
    interpindex = mod(interpindex-2, idealwater.viewnumber) + 1;
    interpindex = index_sort(interpindex);
    % interp rawdata
    dataflow.rawdata = idealwater.main(:, interpindex(:, 1)).*interpalpha(:,1)' +  ...
                       idealwater.main(:, interpindex(:, 2)).*interpalpha(:,2)';
    % interp indexrange
    indexrange = reshape(idealwater.indexrange, 2*idealwater.Nslice, idealwater.viewnumber);
    dataflow.rawhead.index_range = round(indexrange(:, interpindex(:, 1)).*interpalpha(:,1)' + ...
                                   indexrange(:, interpindex(:, 2)).*interpalpha(:,2)');
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
