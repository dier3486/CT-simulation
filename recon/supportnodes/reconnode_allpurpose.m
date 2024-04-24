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
        status.jobdone = true;
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