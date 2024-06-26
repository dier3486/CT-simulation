function [dataflow, prmflow, status] = reconnode_AirAxialcorr(dataflow, prmflow, status)
% recon node, air correction
% [dataflow, prmflow, status] = reconnode_AirAxialcorr(dataflow, prmflow, status);

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

% multi shots
% Nshot = prmflow.raw.Nshot;
if pipeline_onoff
    statuscurr = status.pipepool.(nodename);
    % input shots
    Shot_Number = dataflow.pipepool.(nodename).rawhead.Shot_Number(statuscurr.ReadPoint : statuscurr.WritePoint-1);
    % to check is any data to be output
    % nextnode
    nextnode = status.pipeline.(nodename).nextnode;
    % private buffer ReadPoint
    bufferReadPoint = dataflow.buffer.(nodename).ReadPoint;
    % datasize to write
    datasize_towrite = dataflow.buffer.(nodename).AvailPoint - bufferReadPoint + 1;
    if datasize_towrite > 0
        % we have something to write
        buffersize = dataflow.buffer.(nodename).buffersize;
        % I know buffersize==Nviewprot.
        isfilledaxial = dataflow.buffer.(nodename).AvailPoint == buffersize;
        if isempty(Shot_Number) || isfilledaxial
            % to write to next node
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
            % to recycle private buffer while this
            if isfilledaxial && writenum==datasize_towrite
                % clear private buffer
                dataflow.buffer.(nodename).outpool = ...
                    bufferrecycle(dataflow.buffer.(nodename).outpool, buffersize+1);
                % initial the private Points
                dataflow.buffer.(nodename).ReadPoint = 1;
                dataflow.buffer.(nodename).WritePoint = buffersize - prmflow.raw.air.blockwindow + 1;
                dataflow.buffer.(nodename).AvailPoint = 0;
                % continue as normal (next shot)
                status.jobdone = 1;
                % I konw if the Shot_Number is empty the function will edn soon, if not the next shot will be done.
            else
                % move the private Points
                dataflow.buffer.(nodename).ReadPoint = bufferReadPoint + writenum;
                % job done and keep waking
                status.jobdone = 2;
                status.errorcode = 0;
                status.errormsg = [];
                return;
            end
        else
            % The current shot is not done yet, and there are new data coming and some old data readied in buffer in waiting 
            % to output to next node.
            % It is a dog day but going on as normal.
            1;
        end
    end
else
    Shot_Number = dataflow.rawhead.Shot_Number;
end

shots = unique(Shot_Number);
for ishot = 1:length(shots)
    % shot index
    shotindex = shots(ishot);
    % pipeline consol
    dataflow_redirect = struct();
    if pipeline_onoff
        statuscurr = status.pipepool.(nodename);
        % datasize to read
        datasize_toread = sum(Shot_Number == shotindex);
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
    elseif ishot==1
        % in non-pineline mode do aircorr for all shots once
        dataflow = aircorrwithoutref(dataflow, prmflow, aircorr, sliceindependent);
    end

    if pipeline_onoff
        % back up the refblock
        refblock_old = dataflow.buffer.(nodename).outpool.rawhead.refblock;
        % I know datasize_toread == size(dataflow_redirect.rawdata, 2);
        % copy the dataflow_redirect to (the end of) buffer.outpool
        buffersize = dataflow.buffer.(nodename).poolsize;
        % I know buffersize == Nviewprot.
        redReadPoint = 1;
        bufferWritePoint = dataflow.buffer.(nodename).WritePoint;
        datasize_tobuffer = min(buffersize-bufferWritePoint+1, datasize_toread);
        [dataflow.buffer.(nodename).outpool, ~] = pooldatacopy(dataflow_redirect, dataflow.buffer.(nodename).outpool, ...
            redReadPoint, bufferWritePoint, datasize_tobuffer, [], true);
        bufferWritePoint = bufferWritePoint + datasize_tobuffer;
        if datasize_toread > datasize_tobuffer
            bufferWritePoint = 1;
            redReadPoint = datasize_tobuffer + 1;
            [dataflow.buffer.(nodename).outpool, ~] = pooldatacopy(dataflow_redirect, dataflow.buffer.(nodename).outpool, ...
                redReadPoint, bufferWritePoint, datasize_toread - datasize_tobuffer, [], true);
            bufferWritePoint = datasize_toread - datasize_tobuffer + 1;
        end
        % record the number in writing
        dataflow.buffer.(nodename).Nrenew = datasize_toread;
        % do not over write the refblock
        dataflow.buffer.(nodename).outpool.rawhead.refblock = refblock_old;
    elseif ishot==1
        % ini for normal mode
        %     startrefviewindex = 1;
        %     reflast = cell(1, Nfocal);
        %     dataflow.rawhead.refblock = false(size(altref.referr));
        dataflow.rawhead.refblock = [];
    end

    % air correction step2 (reference correction)
    if pipeline_onoff
        [dataflow.buffer.(nodename).outpool, reflast, Nrenew] = ...
            AirrefcorrAxialKernelfuntion(dataflow.buffer.(nodename).outpool, prmflow, status, shotindex, ...
            dataflow.buffer.(nodename));
    else
        dataflow = AirrefcorrAxialKernelfuntion(dataflow, prmflow, status, shotindex);
    end

    % pipeline consol
    if pipeline_onoff
        % backup the reflast in private buffer
        dataflow.buffer.(nodename).reflast = reflast;
        % move the private buffer Points
        dataflow.buffer.(nodename).AvailPoint = dataflow.buffer.(nodename).AvailPoint + Nrenew;
        dataflow.buffer.(nodename).AvailViewindex = dataflow.buffer.(nodename).AvailViewindex + Nrenew;
        dataflow.buffer.(nodename).WritePoint = bufferWritePoint;

        % nextnode
        nextnode = status.pipeline.(nodename).nextnode;
        % private buffer ReadPoint
        bufferReadPoint = dataflow.buffer.(nodename).ReadPoint;
        % datasize to write
        datasize_towrite = max(0, dataflow.buffer.(nodename).AvailPoint - bufferReadPoint + 1);
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
                dataflow.pipepool.(nextnode), bufferReadPoint, statusnext.WritePoint, writenum);
            % move next pool's write point
            status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
        end
        
        % is this?
        islastshot = ishot == length(shots);
        isfilledaxial = dataflow.buffer.(nodename).AvailPoint == buffersize;

        % to recycle private buffer while this
        if isfilledaxial && writenum==datasize_towrite
            % clear private buffer
            dataflow.buffer.(nodename).outpool = ...
                bufferrecycle(dataflow.buffer.(nodename).outpool, buffersize+1);
            % initial the private Points
            dataflow.buffer.(nodename).ReadPoint = 1;
            dataflow.buffer.(nodename).WritePoint = buffersize - prmflow.raw.air.blockwindow + 1;
            dataflow.buffer.(nodename).AvailPoint = 0;
        else
            % move the private Points
            dataflow.buffer.(nodename).ReadPoint = bufferReadPoint + writenum;
        end
        
        % job done or not
        if datasize_towrite == 0
            if islastshot && ~isfilledaxial
                % did nothing, pass
                status.jobdone = 3;
            else
                % error!
                status.jobdone = false;
                if isfilledaxial
                    status.errorcode = 123;
                    status.errormsg = 'The air correction is not correctly initialized or recycled.';
                else
                    status.errorcode = 113;
                    status.errormsg = 'Not enough views to fill up an axial.';
                end
                return;
            end
        elseif writenum < datasize_towrite
            % done and keep waking
            status.jobdone = 2;
            if isfilledaxial
                % pass the left shots
                break;
            end
        else
            % writenum == datasize_towrite
            if islastshot
                % normally done
                status.jobdone = 1;
            elseif ~isfilledaxial
                % error!
                status.jobdone = false;
                status.errorcode = 112;
                status.errormsg = 'Not enough views to fill up an axial.';
                return;
            end
        end
    else
        status.jobdone = true;
    end

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
    bufferpool.rawhead.(headfields{ii})(:, movesize+1:end) = 0;
end

end
