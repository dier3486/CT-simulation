function [dataflow, prmflow, status] = reconnode_AirHelicalcorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_AirHelicalcorr(dataflow, prmflow, status);

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

nodename = status.nodename;
% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% parameters set in pipe
airprm = prmflow.pipe.(nodename);
% slice independent refernece?
if isfield(airprm, 'sliceindependent')
    sliceindependent = airprm.sliceindependent;
else
    sliceindependent = false;
end

% calibration table
aircorr = prmflow.corrtable.(status.nodename);

% pipeline consol
dataflow_redirect = struct();
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    % datasize to read
    datasize_toread = max(0, statuscurr.WritePoint - statuscurr.ReadPoint);
    % copy pool to dataflow_redirect
    [dataflow_redirect, ~] = pooldatacopy(dataflow.pipepool.(nodename), ...
        dataflow_redirect, statuscurr.ReadPoint, 1, datasize_toread, {}, 1);
    % suppose the private buffer size is infinite.
    % move current pool's read point
    status.pipepool.(nodename).ReadPoint = status.pipepool.(nodename).ReadPoint + datasize_toread; 
end

% air correction step1 ("raw-airtable" and "raw ref") 
if pipeline_onoff
    dataflow_redirect = aircorrwithoutref(dataflow_redirect, prmflow, aircorr, sliceindependent);
else
    dataflow = aircorrwithoutref(dataflow, prmflow, aircorr, sliceindependent);
end

if pipeline_onoff
    % back up the refblock
    refblock_old = dataflow.buffer.(nodename).outpool.rawhead.refblock;
    % I know datasize_toread == size(dataflow_redirect.rawdata, 2);
    % copy the dataflow_redirect to (the end of) buffer.outpool
    WritePoint = dataflow.buffer.(nodename).WritePoint;
    [dataflow.buffer.(nodename).outpool, ~] = pooldatacopy(dataflow_redirect, dataflow.buffer.(nodename).outpool, ...
        1, WritePoint, datasize_toread, [], true);
    % record the number in writing
    dataflow.buffer.(nodename).Nrenew = datasize_toread;
    % do not over write the refblock
    dataflow.buffer.(nodename).outpool.rawhead.refblock = refblock_old;

else
    % ini for normal mode
%     startrefviewindex = 1;
%     reflast = cell(1, Nfocal);
%     dataflow.rawhead.refblock = false(size(altref.referr));
    dataflow.rawhead.refblock = [];
end

% air correction step2 (reference correction)
if pipeline_onoff
    [dataflow.buffer.(nodename).outpool, reflast, Nrenew] = ...
        AirrefcorrHelicalKernelfuntion(dataflow.buffer.(nodename).outpool, prmflow, status, dataflow.buffer.(nodename));
else
    dataflow = AirrefcorrHelicalKernelfuntion(dataflow, prmflow, status);
end

% pipeline consol
if pipeline_onoff
    % backup the reflast in private buffer
    dataflow.buffer.(nodename).reflast = reflast;
    % move the private buffer Points
    dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + Nrenew;
    dataflow.buffer.(nodename).AvailViewindex = dataflow.buffer.(nodename).AvailViewindex + Nrenew;
    dataflow.buffer.(nodename).WritePoint = dataflow.buffer.(nodename).WritePoint + datasize_toread;

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
            dataflow.pipepool.(nextnode), 1, statusnext.WritePoint, writenum);
        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
    end

    % recycle private buffer (must do)
    dataflow.buffer.(nodename).outpool = ...
        bufferrecycle(dataflow.buffer.(nodename).outpool, privReadPoint+writenum);
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
status.errorcode = 0;
status.errormsg = [];
end


function bufferpool = bufferrecycle(bufferpool, recyclePoint)
% buffer recycle

movesize = size(bufferpool.rawdata, 2) - recyclePoint + 1;
bufferpool.rawdata(:, 1:movesize) = bufferpool.rawdata(:, recyclePoint : end);
bufferpool.rawdata(:, movesize+1:end) = 0;

headfields = fieldnames(bufferpool.rawhead);
for ii = 1:length(headfields)
    movesize = size(bufferpool.rawhead.(headfields{ii}), 2) - recyclePoint + 1;
    bufferpool.rawhead.(headfields{ii})(:, 1:movesize) = bufferpool.rawhead.(headfields{ii})(:, recyclePoint : end);
%     bufferpool.rawhead.(headfields{ii})(:, movesize+1:end) = 0;
end

end
