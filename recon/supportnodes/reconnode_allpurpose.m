function [dataflow, prmflow, status] = reconnode_allpurpose(dataflow, prmflow, status)
% support node, all perpose glue
% [dataflow, prmflow, status] = reconnode_allpurpose(dataflow, prmflow, status);
% It is the master node on pipeline used to simulate the function nodes

% prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_allpurposeprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

nodefunct = nodeprm.funct;
nodelevel = nodeprm.level;
echo_onoff = nodeprm.pipeinfoecho;

% if to check data in output pool in start
if isfield(nodeprm, 'checkoutputdatainstart')
    checkoutput_flag = nodeprm.checkoutputdatainstart;
else
    checkoutput_flag = true;
end

% % to stuck in current shot
% if isfield(nodeprm, 'shotstuck')
%     shotstuck_flag = nodeprm.shotstuck;
% else
%     shotstuck_flag = true;
% end

% view common in inputting
if isfield(nodeprm, 'inputviewcommon')
    inputviewcommon = nodeprm.inputviewcommon;
else
    inputviewcommon = 1;
end

% view common in outputting
if isfield(nodeprm, 'outputviewcommon')
    outputviewcommon = nodeprm.outputviewcommon;
else
    outputviewcommon = 1;
end

% % view rely
% if isfield(nodeprm, 'viewrely')
%     viewrely = nodeprm.viewrely;
% else
%     viewrely = [0 0];
% end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;


if pipeline_onoff
    % next node
    nextnode = status.pipeline.(nodename).nextnode;
    if isempty(nextnode)
        nextnode = 'NULL';
    end 
    % echo
    if echo_onoff
        if ~isempty(nodefunct)
            fprintf(' (%s)', nodefunct);
        end
        fprintf('\n        ');
    end
end

% prio
if pipeline_onoff
    % prio step #1
    % try to output output pool to next pool
    if checkoutput_flag && nodelevel>=1
        % output pool is dataflow.buffer.(nodename).outputpool
        % next pool is dataflow.pipepool.(nextnode)
        1;
        [dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, status.pipepool.(nextnode), ...
            dataflow.pipepool.(nextnode), priojobdone1] = ...
            clearoutputarray(dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
            status.pipepool.(nextnode), dataflow.pipepool.(nextnode), outputviewcommon, echo_onoff);
        if priojobdone1 >= 2
            % stuck due to next pool is filled
            status.jobdone = priojobdone1;
            return;
        end
    else
        priojobdone1 = 0;
    end
    
    % prio step #2
    % to get data from current pool
    switch nodelevel
        case 0
            % level 0, copy current pool to next pool
            % current pool is dataflow.pipepool.(nodename)
            % next pool is dataflow.pipepool.(nextnode)
            1;
            [status.pipepool.(nodename), status.pipepool.(nextnode), dataflow.pipepool.(nextnode), ...
                inputRinfo] = pooldataoutput2(status.pipepool.(nodename), dataflow.pipepool.(nodename), ...
                status.pipepool.(nextnode), dataflow.pipepool.(nextnode), ...
                inputviewcommon);

        case 1
            1;
            % level 1, copy current pool to output pool
            % current pool is dataflow.pipepool.(nodename)
            % output pool is dataflow.buffer.(nodename).outputpool
            [status.pipepool.(nodename), dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                inputRinfo] = pooldataoutput2(status.pipepool.(nodename), dataflow.pipepool.(nodename), ...
                dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                inputviewcommon);
            1;
            
        case {2, 3}
            % level 2, copy current pool to input pool
            % current pool is dataflow.pipepool.(nodename)
            % input pool is dataflow.buffer.(nodename).inputpool
            [status.pipepool.(nodename), dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).inputpool.data, ...
                inputRinfo] = pooldataoutput2(status.pipepool.(nodename), dataflow.pipepool.(nodename), ...
                dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).inputpool.data, ...
                inputviewcommon);
        otherwise
            % ??
            inputRinfo.writenum = 0;
    end
    if echo_onoff
        fprintf('read from previous %d views, ', inputRinfo.writenum);
    end

    priojobdone2 = inputRinfo.jobflag;
    
    % Note: do not recycle the current pool!
else
    inputRinfo = struct();
end


% node main
switch nodefunct
    case {'loadrawdata', 'readraw'}
        [dataflow, prmflow, status] = fakenodereadrawdata(dataflow, prmflow, status);
    case 'donothing'
        1;
    otherwise
        % default node function
        1;
        [dataflow, prmflow, status] = allpurposefakenode(dataflow, prmflow, status, inputRinfo);
end

% post
if pipeline_onoff
    % post step #1
    % to output data to next pool
    

    switch nodelevel
        case 0
            % level 0
            % move nextpool's AvailNumber
            if isfield(status.pipepool.(nextnode), 'AvailNumber')
                status.pipepool.(nextnode).AvailNumber = status.pipepool.(nextnode).AvailNumber + inputRinfo.writenum;
            end
            if echo_onoff
                fprintf('write to next %d views, ', inputRinfo.writenum);
            end
            % we have done the data copy
            postjobdone = priojobdone1;
        case 1
            % level 1, output output pool to next pool and recycle the output pool
            % output pool is dataflow.buffer.(nodename).outputpool
            % next pool is dataflow.pipepool.(nextnode)
            [dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                status.pipepool.(nextnode), dataflow.pipepool.(nextnode), postjobdone] = ...
                clearoutputarray(dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                status.pipepool.(nextnode), dataflow.pipepool.(nextnode), outputviewcommon, echo_onoff);
        case {2, 3}
            % level 2
            % step #1.1 recycle the input pool
            [dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).inputpool.data] = ...
                poolrecycle(dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).inputpool.data);
            
            1;
            % step #1.2 output output pool to next pool and recycle the output pool
            [dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                status.pipepool.(nextnode), dataflow.pipepool.(nextnode), postjobdone] = ...
                clearoutputarray(dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data, ...
                status.pipepool.(nextnode), dataflow.pipepool.(nextnode), outputviewcommon, echo_onoff);
            1;
        otherwise
    end
    
    % post step #2
    % to return the status
    [status.jobdone, errorcode] = jobdonetable(priojobdone1, priojobdone2, status.jobdone, postjobdone);
    if status.jobdone==0 && ~strcmp(errorcode, '0')
        status.errorcode = ['3' errorcode];
        status.errormsg = ['Pipe-line error event ' errorcode];
    end
    
end

end


function [dataflow, prmflow, status] = fakereadraw(dataflow, prmflow, status)

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% shots to read
shotnum = prmflow.raw.Nshot;
startshot = prmflow.raw.startshot;
viewpershot = prmflow.raw.viewpershot;
if isfield(dataflow, 'buffer') && dataflow.buffer.(nodename).datablock_onoff
    % to read data by blocks
    datablocksize = dataflow.buffer.(nodename).datablocksize;
    iblock = dataflow.buffer.(nodename).iblock;
    if iblock > length(datablocksize)
        % I think it is a mistake in pipeline
        % pass
        status.jobdone = 3;
        return
    end
    startview = sum(datablocksize(1:iblock-1)) + (startshot-1)*viewpershot + 1;
    viewnum = datablocksize(iblock);
else
    startview = (startshot-1)*viewpershot + 1;
    viewnum = shotnum * viewpershot;
end

% fake rawdata
if ~isfield(dataflow, 'rawdata')
    dataflow.rawdata = [];
end
dataflow.rawdata = [dataflow.rawdata, zeros(1, viewnum)];
if ~isfield(dataflow, 'rawhead')
    dataflow.rawhead = struct();
    dataflow.rawhead.Reading_Number = [];
    dataflow.rawhead.Shot_Number = [];
end
currReading_Number = startview : startview+viewnum-1;
currShotnumber = ceil(currReading_Number / viewpershot);
dataflow.rawhead.Reading_Number = [dataflow.rawhead.Reading_Number currReading_Number];
dataflow.rawhead.Shot_Number = [dataflow.rawhead.Shot_Number currShotnumber];

% datablock
if isfield(dataflow, 'buffer') && dataflow.buffer.(nodename).datablock_onoff
    % contol the looping of data blocks
    iblock = dataflow.buffer.(nodename).iblock;
    Nblock = dataflow.buffer.(nodename).Nblock;
    if iblock < Nblock
        status.jobdone = 2;
    else
%         status.pipeline.(nodename).sleeping = true;
        status.jobdone = 1;
    end
    % iblock++
    dataflow.buffer.(nodename).iblock = iblock + 1;
    % WritePoint for the 'hidden' node
    if isfield(dataflow, 'buffer')
        dataflow.buffer.(nodename).pipepool.WritePoint = dataflow.buffer.(nodename).pipepool.WritePoint + viewnum;
    end
else
    status.jobdone = true;
end

if pipeline_onoff && ~dataflow.buffer.(nodename).datablock_onoff
    % What? Shouldn't we forbid this behavior?
    dataflow.buffer.(nodename).pipepool.WritePoint = dataflow.buffer.(nodename).pipepool.WritePoint + viewnum;
end

if pipeline_onoff
% pipeline hidden node
if pipeline_onoff
    nextnode = status.pipeline.(nodename).nextnode;
    % the current pool of the 'hidden' node is not in status.pipepool, but in private buffer
    currpool = dataflow.buffer.(nodename).pipepool;
    if ~isempty(nextnode)
        statusnext = status.pipepool.(nextnode);
        % data size to be written in the next pool
        writenum = min(viewnum, min(statusnext.WriteEnd, statusnext.poolsize) - statusnext.WritePoint + 1);
        % copy rawdata and rawhead to next pool
        1;
        [dataflow.pipepool.(nextnode), writenum] = pooldatacopy(dataflow, dataflow.pipepool.(nextnode), ...
            currpool.ReadPoint, statusnext.WritePoint, writenum, {'rawdata', 'rawhead'}, true);

        % move next pool's write point
        status.pipepool.(nextnode).WritePoint = status.pipepool.(nextnode).WritePoint + writenum;
        status.pipepool.(nextnode).AvailNumber = status.pipepool.(nextnode).AvailNumber + writenum;
        % Maybe we shall let pipeline consol doing that.
        % move current pool's read point
        dataflow.buffer.(nodename).pipepool.ReadPoint = dataflow.buffer.(nodename).pipepool.ReadPoint + writenum;
        % recycle
        currpool = dataflow.buffer.(nodename).pipepool;
        % datalength
        datalength = currpool.WritePoint - currpool.ReadPoint;
        switch currpool.recylestrategy
            case 0
                recycle_onoff = false;
            case 1
                % recycle
                recycle_onoff = true;
            case 2
                % recycle?
                recycle_onoff = currpool.WritePoint > currpool.warningstage;
            otherwise
                status.jobdone = false;
                status.errorcode = 901;
                status.errormsg = sprintf('readrawdata''s hidden node error, unkown recylestrategy %d!', currpool.recylestrategy);
                return
        end
        if recycle_onoff && datalength>0
            dataflow.rawdata(:, 1:datalength) = ...
                dataflow.rawdata(:, currpool.ReadPoint : currpool.WritePoint-1);
            headfields = fieldnames(dataflow.rawhead);
            for ii = 1:length(headfields)
                dataflow.rawhead.(headfields{ii})(:, 1:datalength) = ...
                    dataflow.rawhead.(headfields{ii})(:, currpool.ReadPoint : currpool.WritePoint-1);
            end
            dataflow.buffer.(nodename).pipepool.ReadPoint = 1;
            dataflow.buffer.(nodename).pipepool.WritePoint = datalength + 1;
        end
        
    end 
end
end

end