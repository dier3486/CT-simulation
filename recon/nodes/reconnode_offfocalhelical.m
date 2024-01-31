function [dataflow, prmflow, status] = reconnode_offfocalhelical(dataflow, prmflow, status)
% recon node, off-focal correction for helical
% [dataflow, prmflow, status] = reconnode_offfocalhelical(dataflow, prmflow, status);

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

% nodename
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipeline consol
% read the data from input pool
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    % datasize
    datasize_toread = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % check the buffer left
    outpoolsize = dataflow.buffer.(nodename).outpoolsize - dataflow.buffer.(nodename).WritePoint + 1;
    datasize_toread = min(datasize_toread, outpoolsize);

    % copy pool to buffer.outpool
    [dataflow.buffer.(nodename).outpool, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow.buffer.(nodename).outpool, statuscurr.ReadPoint, dataflow.buffer.(nodename).WritePoint, ...
        datasize_toread, {}, true);

    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + datasize_toread; 
    % move the buffer.outpool's write point
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint + datasize_toread;
    % Note: the buffer.(nodename).outpool is a private buffer, not the next pool. The data in buffer.(nodename).outpool will be
    % copied to the next pool after the correcion works done.
end

% step1 exp
if pipeline_onoff
    % view index to do 
    index = (1:datasize_toread) - datasize_toread + dataflow.buffer.(nodename).WritePoint - 1;
    % the dataflow is redirected to dataflow.buffer.(nodename).outpool
    dataflow.buffer.(nodename).outpool.rawdata(:, index) = 2.^(-dataflow.buffer.(nodename).outpool.rawdata(:, index));
else
    dataflow.rawdata = 2.^(-dataflow.rawdata);
end

% step2, resample to off-focal space and convolution
if pipeline_onoff
    % the dataflow is redirected to dataflow.buffer.(nodename).outpool
    offraw = offfocalmeasurehelicalKernelfunction(dataflow.buffer.(nodename).outpool, prmflow, status, ...
        dataflow.buffer.(nodename));
else
    offraw = offfocalmeasurehelicalKernelfunction(dataflow, prmflow, status);
end

% pipeline consol
% copy the offraw to buffer.offfocalspace
if pipeline_onoff
    % data size to be written in offfocalspace
    writenum = size(offraw.rawdata, 2);
    % suppose the size of buffer.offfocalspace is Inf. 
    % (If not, plz modify the codes in offfocalmeasurehelicalKernelfunction, in working...)

    % copy offraw.rawdata to offfocalspace.rawdata
    [dataflow.buffer.(nodename).offfocalspace, ~] = pooldatacopy(offraw, dataflow.buffer.(nodename).offfocalspace, ...
        1, dataflow.buffer.(nodename).offWritePoint, writenum, {'rawdata'}, true);

    % move the buffer.offfocalspace's write point
    dataflow.buffer.(nodename).offWritePoint = dataflow.buffer.(nodename).offWritePoint + writenum;
end

% step3, resample back and add to rawdata
if pipeline_onoff
    [dataflow.buffer.(nodename).outpool, Nrenew, offNremove] = offfocalfixhelicalKernelfunction(...
        dataflow.buffer.(nodename).outpool, prmflow, status, dataflow.buffer.(nodename).offfocalspace, ...
        dataflow.buffer.(nodename));
else
    dataflow = offfocalfixhelicalKernelfunction(dataflow, prmflow, status, offraw);
end

% pipeline consol
% copy the buffer.outpool to next pool, recycle the private buffer and move the points
if pipeline_onoff
    % move the outpool's available point and view index
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + Nrenew;
    dataflow.buffer.(nodename).AvailViewindex = dataflow.buffer.(nodename).AvailViewindex + Nrenew;

    % recycle the buffer.offfocalspace
    dataflow.buffer.(nodename).offfocalspace = poolrecycle(dataflow.buffer.(nodename).offfocalspace, offNremove);
    % keep the buffer.offfocalspace's read point = 1,
    % move it's vieindex
    dataflow.buffer.(nodename).offReadViewindex = dataflow.buffer.(nodename).offReadViewindex + offNremove;
    % move the buffer.offfocalspace's write point
    dataflow.buffer.(nodename).offWritePoint = dataflow.buffer.(nodename).offWritePoint - offNremove;

    % I will move these codes to public
    % nextnode
    nextnode = status.pipeline.(nodename).nextnode;
    % private ReadPoint
    privReadPoint = dataflow.buffer.(nodename).ReadPoint;
    % datasize to write
    datasize_towrite = dataflow.buffer.(nodename).AvailPoint - privReadPoint + 1;
    % writenum
    if ~isempty(nextnode)
        statusnext = status.pipepool.(nextnode);
        writenum = min(datasize_towrite, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
    else
        writenum = datasize_towrite;
        % Even when the nextnode is not existing the node will do its work as usual.
    end
    % to next pool
    if ~isempty(nextnode)
        % copy dataflow_redirect to next pool
        [dataflow.pipepool.(nextnode), ~] = pooldatacopy(dataflow.buffer.(nodename).outpool, ...
            dataflow.pipepool.(nextnode), privReadPoint, statusnext.WritePoint, writenum);
        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
    end
    % recycle private buffer (must do)
    removenumber = privReadPoint + writenum - 1;    % max recycle
    dataflow.buffer.(nodename).outpool = ...
        poolrecycle(dataflow.buffer.(nodename).outpool, removenumber);
    % move the private Points
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint - removenumber;
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint - removenumber;

    % job done
    status.jobdone = dailyjobdone(datasize_toread, writenum, datasize_towrite);    
end
1;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end