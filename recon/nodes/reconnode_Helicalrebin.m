function [dataflow, prmflow, status] = reconnode_Helicalrebin(dataflow, prmflow, status)
% recon node, Helical rebin 
% [dataflow, prmflow, status] = reconnode_Helicalrebin(dataflow, prmflow, status);
% Support (X)DFS, gantry tilt, NO QDO

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
    [dataflow, prmflow, status] = reconnode_helicalrebinprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% nodename
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipeline consol
dataflow_redirect = struct();
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
%     nextnode = status.pipeline.(nodename).nextnode;
    % datasize
    datasize_toread = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % only even
    datasize_toread = floor(datasize_toread/2)*2;
    % copy pool to dataflow_redirect
    [dataflow_redirect, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow_redirect, statuscurr.ReadPoint, 1, datasize_toread, {}, 1);
    % suppose the private buffer size is infinite.
    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + datasize_toread; 
end

% step 1 Radial
if pipeline_onoff
    [dataflow_redirect, prmflow] = FanRadialrebinKernelfunction(dataflow_redirect, prmflow, status);
else
    [dataflow, prmflow] = FanRadialrebinKernelfunction(dataflow, prmflow, status);
end

% step2 Azi
if pipeline_onoff
    [dataflow_redirect, prmflow] = HelicalAzirebinKernelfuntion(dataflow_redirect, prmflow, status);
else
    [dataflow, prmflow] = HelicalAzirebinKernelfuntion(dataflow, prmflow, status);
end

% save in private buffer and copy back
if pipeline_onoff
    % rebin prm
    Nview = prmflow.rebin.Nview;
    Nextraview = prmflow.rebin.Nextraview;
    % private WritePoint and redirected ReadPoint
    privWritePoint = dataflow.buffer.(nodename).WritePoint;
    redirectReadPoint = max(0, -dataflow.buffer.(nodename).AvailViewindex) + 1;
    % data cum up
    dataflow.buffer.(nodename).outpool = pooldatacum(dataflow_redirect, dataflow.buffer.(nodename).outpool, ...
        redirectReadPoint, privWritePoint);
    % move the private WritePoint and other
    dataflow.buffer.(nodename).WritePoint = privWritePoint + datasize_toread - redirectReadPoint + 1;
    dataflow.buffer.(nodename).AvailViewindex = dataflow.buffer.(nodename).AvailViewindex + datasize_toread;
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + datasize_toread;
    % if reach the end view
    if dataflow.buffer.(nodename).AvailViewindex >= Nview - Nextraview(1)
        dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + Nextraview(1);
    end

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
    dataflow.buffer.(nodename).outpool = ...
        bufferrecycle(dataflow.buffer.(nodename).outpool, privReadPoint, writenum);
    % move the private Points
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint - writenum;
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint - writenum;

    % job done
    if datasize_towrite == 0 && datasize_toread == 0
        % did nothing, pass
        status.jobdone = 3;
    elseif datasize_towrite == 0
        % read something but write nothing, technically pass
        status.jobdone = 4;
    elseif writenum < datasize_towrite
        % done and keep waking
        status.jobdone = 2;
    else
        % normally done
        status.jobdone = 1;
    end
else
    status.jobdone = true;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

function bufferpool = pooldatacum(dataflow, bufferpool, ReadPoint, WritePoint)
% data cum up is to link two block of data with overlap, in overlap part the two block of data will be sum up
% I will move this function to public.

Nview_write = size(dataflow.rawdata, 2);
writenum = Nview_write - ReadPoint + 1;

if ~isfield(bufferpool, 'rawdata')
    bufferpool.rawdata = dataflow.rawdata(:, ReadPoint : end);
else
    poollength = size(bufferpool.rawdata, 2);
    if poollength >= (WritePoint + writenum - 1)
        bufferpool.rawdata(:, WritePoint : WritePoint + writenum - 1) = ...
            bufferpool.rawdata(:, WritePoint : WritePoint + writenum - 1) + dataflow.rawdata(:, ReadPoint : end);
    else
        Nfill = poollength - WritePoint;
        bufferpool.rawdata(:, WritePoint : end) = bufferpool.rawdata(:, WritePoint : end) + ...
            dataflow.rawdata(:, ReadPoint : ReadPoint + Nfill);
        bufferpool.rawdata = [bufferpool.rawdata  dataflow.rawdata(:, ReadPoint + Nfill + 1 : end)];
    end
end

headfields = fieldnames(dataflow.rawhead);
if ~isfield(bufferpool, 'rawhead')
    bufferpool.rawhead = struct();
    for ii = 1:length(headfields)
        bufferpool.rawhead.(headfields{ii}) = dataflow.rawhead.(headfields{ii})(:, ReadPoint : end);
    end
else
    if poollength >= (WritePoint + writenum - 1)
        for ii = 1:length(headfields)
            bufferpool.rawhead.(headfields{ii})(:, WritePoint : WritePoint + writenum - 1) = ...
                bufferpool.rawhead.(headfields{ii})(:, WritePoint : WritePoint + writenum - 1) + ...
                dataflow.rawhead.(headfields{ii})(:, ReadPoint : end);
        end
    else
        for ii = 1:length(headfields)
            bufferpool.rawhead.(headfields{ii})(:, WritePoint : end) = ...
                bufferpool.rawhead.(headfields{ii})(:, WritePoint : end) + ...
                dataflow.rawhead.(headfields{ii})(:, ReadPoint : ReadPoint + Nfill);
            bufferpool.rawhead.(headfields{ii}) = [bufferpool.rawhead.(headfields{ii})  ...
                dataflow.rawhead.(headfields{ii})(:, ReadPoint + Nfill + 1: end)];
        end
    end
end

end


function bufferpool = bufferrecycle(bufferpool, ReadPoint, writenum)
% to recycle a private buffer
% to be replaced by public function poolrecycle

movesize = size(bufferpool.rawdata, 2) - ReadPoint - writenum + 1;
bufferpool.rawdata(:, 1:movesize) = bufferpool.rawdata(:, ReadPoint + writenum : end);
bufferpool.rawdata(:, movesize+1:end) = 0;

headfields = fieldnames(bufferpool.rawhead);
for ii = 1:length(headfields)
    bufferpool.rawhead.(headfields{ii})(:, 1:movesize) = bufferpool.rawhead.(headfields{ii})(:, ReadPoint + writenum : end);
    bufferpool.rawhead.(headfields{ii})(:, movesize+1:end) = 0;
end

end
