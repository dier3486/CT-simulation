function [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status)
% recon node, filter design and conv
% [dataflow, prmflow, status] = reconnode_filter(dataflow, prmflow, status);

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
    [dataflow, prmflow, status] = reconnode_filterprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% prm
Nslice = prmflow.rebin.Nslice;
nodeprm = prmflow.pipe.(status.nodename);

% filter
basicfilter = prmflow.recon.filter.basicfilter;
Npixel = prmflow.recon.filter.Npixel;
Hlen = prmflow.recon.filter.Hlen;

% fill up
if isfield(nodeprm, 'fillup') && nodeprm.fillup
    % to fill up the data for off-ISOcenter
    if isfield(dataflow.rawhead, 'refblock')
        blkvindex = any(dataflow.rawhead.refblock, 1);
    else
        blkvindex = [];
    end
    
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, []);
    Nview = size(dataflow.rawdata, 3);
    A = zeros(Hlen, Nslice, Nview);
    for ii = 1:Nslice
        [A(:, ii, :), n_left] = translatefillup(squeeze(dataflow.rawdata(:, ii, :)), Hlen, mid_u, blkvindex);
    end
    dataflow.rawdata = reshape(A, Hlen, []);
else
    % fill zero
    dataflow.rawdata = reshape(dataflow.rawdata, Npixel, []);
    dataflow.rawdata(Hlen, :) = 0;
    n_left = 0;
end

% conv
% fft
dataflow.rawdata = fft(dataflow.rawdata);
% timesymmetric
dataflow.rawdata = dataflow.rawdata.*basicfilter;
% ifft
dataflow.rawdata = ifft(dataflow.rawdata, 'symmetric');
% kick filled zero
% dataflow.rawdata(Npixel+1:end,:) = [];
dataflow.rawdata = dataflow.rawdata((1:Npixel)+n_left, :);
% done

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end