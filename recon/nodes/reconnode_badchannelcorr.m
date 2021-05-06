function [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status)
% recon node, interp for bad channels
% [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status);

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

% parameters to use
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;

% calibration table
if isfield(prmflow.corrtable, status.nodename)
    badcorr = prmflow.corrtable.(status.nodename);
    badindex = badcorr.badindex(:);
else
    badindex = [];
end

% input bad channel index
if isfield(prmflow.pipe, status.nodename)
    nodeprm = prmflow.pipe.(status.nodename);
    if isfield(nodeprm, 'badindex')
        if ischar(nodeprm.badindex)
            index_add = cellfun(@str2num, regexp(nodeprm.badindex, '\d+', 'match'));
        else
            index_add = nodeprm.badindex;
        end
        badindex = unique([badindex; index_add(:)]);
    end
end

% find nan
dataflow.rawdata = reshape(dataflow.rawdata, Nps, []);
index_add = find(isnan(sum(dataflow.rawdata, 2)));
badindex = unique([badindex; index_add]);

badmax = 10;
% loop the bad channels
for ii = 1:length(badindex)
    ibad = badindex(ii);
    % to find the left good one
    for jj = 1:badmax
        i_left = ibad-jj;
        if any(i_left==badindex)
            % bad again
            i_left = [];
        elseif mod(i_left, Npixel) == 0
            % touch the edge
            i_left = [];
            a_left = inf;
            break;
        else
            % got
            a_left = jj;
            break;
        end
    end
    
    % to find the right good one
    for jj = 1:badmax
        i_right = ibad+jj;
        if any(i_right==badindex)
            % bad again
            i_right = [];
        elseif mod(i_right, Npixel) == 1
            % touch the edge
            i_right = [];
            a_right = inf;
            break;
        else
            % got
            a_right = jj;
            break;
        end
    end
    % interp weight
    alpha = [a_right, a_left]./(a_right+a_left);
    alpha(isnan(alpha)) = 0;
    % interp
    dataflow.rawdata(ibad, :) = alpha * dataflow.rawdata([i_left i_right], :);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end