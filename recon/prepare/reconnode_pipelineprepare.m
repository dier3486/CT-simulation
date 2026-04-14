function [dataflow, prmflow, status] = reconnode_pipelineprepare(dataflow, prmflow, status)
% pipeline prepare node
% [dataflow, prmflow, status] = reconnode_pipelineprepare(dataflow, prmflow, status);
% Plz run this node after loading calibration tables

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

% echo onoff
echo_onoff = status.debug.echo_onoff;

% backup status.nodename
orignodename = status.nodename;

% % prepare of read rawdata
% prmflow = readrawdataprepare(prmflow, status);
% status.pipeline.loadrawdata.prepared = true;

% prepare of the pipeline nodes
pipenodes = fieldnames(status.pipeline);

flag_prepared = 0;
% to record the error message from the prepare fucntions 
preperrormsg = '';
preperrorcode = 0;
% loop the nodes
for ii = 1:length(pipenodes)
    % spell the nodes prepare function
    nodename_slip = regexp(pipenodes{ii}, '_', 'split');
    preparename = ['reconnode_' lower(nodename_slip{1}) 'prepare'];
    % Do not name a pipe line node in 'pipeline' ;)  
    status.nodename = pipenodes{ii};
    status.currentjob.prepared = false;  
    % that flag is used to close the defaultpipelineprepare, while a node's
    % prepare function set it to true the default-prepare will be skiped.
    % to call the node's prepare function
    if any(exist(preparename) == [2 5 6]) % can run
        if echo_onoff
            if flag_prepared > 0
                fprintf(', ');
            else
                fprintf(' [');
            end
            flag_prepared = flag_prepared + 1;
            if mod(flag_prepared, 6) == 0
                fprintf('\n\t\t');
            end
            fprintf(pipenodes{ii});
        end
        preparefun = str2func(preparename);
        [dataflow, prmflow, status] = preparefun(dataflow, prmflow, status);
        if status.jobdone == 0
            return
        elseif status.errorcode ~=0
            % record them
            preperrorcode = status.errorcode;
            preperrormsg = [preperrormsg status.errormsg ' '];
            status.errorcode = 0;
            status.errormsg = '';
        end
        % set the flag prepared in status.pipeline to true
        status.pipeline.(status.nodename).prepared = true;
    end
    % default prepares of pipepool
    if ~status.currentjob.prepared
        [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status);
    end
end
if echo_onoff && flag_prepared > 0
    fprintf(']');
end

% auto set carry nodes, only for pipepool(1) no branches
nodename = pipenodes{length(pipenodes)};
while ~strcmpi(nodename, 'NULL')
    % prevnode
    prevnode = status.pipeline.(nodename).prevnode;
    if strcmpi(prevnode, 'NULL')
        % done while meating a NULL node
        break;
        % continue if in these condition
    elseif ~status.pipeline.(prevnode).pipeline_onoff
        nodename = prevnode;
        continue;
    elseif ~status.pipeline.(nodename).pipeline_onoff
        nodename = prevnode;
        continue;
    end
    % check the carry of current node
    if prmflow.pipe.(nodename).pipeline.iscarried && ~dataflow.pipepool.(nodename)(1).iscarried
        prmflow.pipe.(nodename).pipeline.iscarried = false;
    end

    % is (the previous node) carry
    if prmflow.pipe.(prevnode).pipeline.iscarried
        % carried
        dataflow.pipepool.(prevnode)(1).iscarried = true;
        % align these,
        prmflow.pipe.(prevnode).pipeline.currcirculte = prmflow.pipe.(nodename).pipeline.currcirculte;
        if prmflow.pipe.(prevnode).pipeline.inputminlimit < prmflow.pipe.(nodename).pipeline.inputminlimit
            prmflow.pipe.(prevnode).pipeline.inputminlimit = prmflow.pipe.(nodename).pipeline.inputminlimit;
        end
        dataflow.pipepool.(prevnode)(1).circulatemode = dataflow.pipepool.(nodename)(1).circulatemode;
        setpoolsize = max(dataflow.pipepool.(prevnode)(1).poolsize, dataflow.pipepool.(nodename)(1).poolsize);
        dataflow.pipepool.(prevnode)(1).poolsize = setpoolsize;
        dataflow.pipepool.(nodename)(1).poolsize = setpoolsize;
        dataflow.pipepool.(prevnode)(1).data = struct();  % not []
        if ~strcmp(prmflow.pipe.(prevnode).pipeline.objecttype, prmflow.pipe.(nodename).pipeline.objecttype)
            preperrorcode = 231;
            errormsg = ...
                sprintf('Node %s and %s, can not carry the nodes in different objecttype.' ...
                , prevnode, nodename);
            preperrormsg = [preperrormsg errormsg ' '];
%         elseif prmflow.pipe.(prevnode).pipeline.kernellevel ~= 0
%             preperrorcode = 231;
%             errormsg = ...
%                 sprintf('Node %s, can not carry the node with kernellevel>0.' ...
%                 , prevnode);
%             preperrormsg = [preperrormsg errormsg ' '];
        end
        % carry pool
        if prmflow.pipe.(nodename).pipeline.iscarried
            dataflow.pipepool.(prevnode)(1).carrynode = dataflow.pipepool.(nodename)(1).carrynode;
%             status.pipeline.(prevnode).carrynode = status.pipeline.(nodename).carrynode;
        else
            dataflow.pipepool.(prevnode)(1).carrynode = nodename;
%             status.pipeline.(prevnode).carrynode = nodename;
        end
        % I know the default carryindex is 1, need not change.       
    end
    nodename = prevnode;
end

% add carried pools to carriages
for ii = 1:length(pipenodes)
    nodename = pipenodes{ii};
    if strcmpi(nodename, 'NULL') || ~status.pipeline.(nodename).pipeline_onoff
        continue;
    end
    for jj = 1:length(dataflow.pipepool.(nodename))
        if dataflow.pipepool.(nodename)(jj).iscarried
            carrynode = dataflow.pipepool.(nodename)(jj).carrynode;
            carryindex = dataflow.pipepool.(nodename)(jj).carryindex;
            % check if the carrypool is final
            if dataflow.pipepool.(carrynode)(carryindex).iscarried
                % to find out the final carrypool
                [carrynode, carryindex] = carryrecurse(dataflow.pipepool, carrynode, carryindex);
                dataflow.pipepool.(nodename)(jj).carrynode = carrynode;
                dataflow.pipepool.(nodename)(jj).carryindex = carryindex;
            end
            if ~any(strcmp(dataflow.pipepool.(carrynode)(carryindex).carriages, nodename))
                dataflow.pipepool.(carrynode)(carryindex).carriages = ...
                    [dataflow.pipepool.(carrynode)(carryindex).carriages nodename];
                1;
            end
        end
    end
end

% double check the prepared nodes (to return error), only for debug
for ii = 1:length(pipenodes)
    nodename = pipenodes{ii};
    if strcmpi(nodename, 'NULL') || ~status.pipeline.(nodename).pipeline_onoff
        continue;
    end
    % check an illegal configure of the kernellevel 0
    if prmflow.pipe.(nodename).pipeline.kernellevel==0 && any(prmflow.pipe.(nodename).pipeline.viewrescale ~= 1) ...
            && prmflow.pipe.(nodename).pipeline.viewrescale(1) ~= prmflow.pipe.(nodename).pipeline.viewrescale(2)
        prmflow.pipe.(nodename).pipeline.kernellevel = 1;
        preperrorcode = 232;
        errormsg = sprintf('Node %s, the kernel level can not be 0 while viewrescaling! Which is forcedly set to 1.'...
            , nodename);
        preperrormsg = [preperrormsg errormsg ' '];
    end
    % nextnode
    nextnode = status.pipeline.(nodename).nextnode;
    prevnode = status.pipeline.(nodename).prevnode;
    % word nodetype
    circulateS = 'HA';
    levelS = '01';
    strategyS = 'NGS';
    if strcmpi(nextnode, 'NULL') || ~status.pipeline.(nextnode).pipeline_onoff
        prmflow.pipe.(nodename).pipeline.nodetype = 'NULL';
        % I know the last node's nextnode is NULL, and a node laying on pipeline_onoff==false will be skipped in
        % dataflow.pipepool
    else
        prmflow.pipe.(nodename).pipeline.nodetype = [circulateS(dataflow.pipepool.(nodename)(1).circulatemode+1) '-' ...
            circulateS(dataflow.pipepool.(nextnode)(1).circulatemode+1) '.' ...
            levelS(prmflow.pipe.(nodename).pipeline.kernellevel+1) '.' ...
            strategyS(prmflow.pipe.(nodename).pipeline.relystrategy+1)];
    end
    % check illegal nodetype
    if strcmp(prmflow.pipe.(nodename).pipeline.nodetype, 'H-A.1.S')
        preperrorcode = 232;
        errormsg = ...
            sprintf('Node %s: illegal nodetype H-A.1.S.' ...
            , nodename);
        preperrormsg = [preperrormsg errormsg ' '];
    end
    % % check illegal view extra
    % if ~strcmpi(nextnode, 'NULL')
    %     if ~strcmpi(prmflow.pipe.(nodename).pipeline.nodetype, 'NULL') && status.pipeline.(nextnode).pipeline_onoff
    %         if any(prmflow.pipe.(nodename).pipeline.viewextra ~= 0) && (dataflow.pipepool.(nodename)(1).circulatemode || ...
    %                 dataflow.pipepool.(nextnode)(1).circulatemode)
    %             prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
    %             preperrorcode = 232;
    %             errormsg = ...
    %                 sprintf('Node %s, can not employ the view-extra in circulate-mode! Which is forcedly set to [0 0].' ...
    %                 , nodename);
    %             preperrormsg = [preperrormsg errormsg ' '];
    %         end
    %     end
    % end
    % check carry node
    %         if prmflow.pipe.(nodename).pipeline.iscarried && prmflow.pipe.(nodename).pipeline.kernellevel>0
    %             preperrorcode = 232;
    %             errormsg = ...
    %                     sprintf('Node %s, only a level 0 node can be carried.' ...
    %                     , nodename);
    %                 preperrormsg = [preperrormsg errormsg ' '];
    %                 % yes, we can
    %         end
    % check poolsize
    if ~strcmpi(prevnode, 'NULL')
        pooltouse = prmflow.pipe.(nodename).pipeline.inputminlimit + sum(prmflow.pipe.(prevnode).pipeline.viewrely_out);
    else
        pooltouse = prmflow.pipe.(nodename).pipeline.inputminlimit;
    end
    if ~dataflow.pipepool.(nodename)(1).circulatemode && dataflow.pipepool.(nodename)(1).poolsize < pooltouse
        preperrorcode = 234;
        errormsg = ...
            sprintf('Node %s, the poolsize is too small in probable stucking risk.' ...
            , nodename);
        preperrormsg = [preperrormsg errormsg ' '];
    end
    % Note: not all stucking risks can be found by that.
    % check GPU
%     objecttype = prmflow.pipe.(nodename).pipeline.objecttype;
%     if ~dataflow.pipepool.(nodename)(1).iscarried
%         if (prmflow.pipe.(nodename).pipeline.currgpuDevice > 0) ~= isgpuarray(dataflow.pipepool.(nodename)(1).data.(objecttype))
%             preperrorcode = 234;
%             errormsg = ...
%             sprintf('Node %s, the pool.data is not been correctly claimed by GPU on/off.' ...
%             , nodename);
%         preperrormsg = [preperrormsg errormsg ' '];
%         end
%     end

end

% other things 
for ii = 1:length(pipenodes)
    nodename = pipenodes{ii};
    if strcmpi(nodename, 'NULL')
        continue;
    end
    % record the GPU device in status (maybe useless)
    if prmflow.pipe.(nodename).pipeline.GPUonoff
        if prmflow.pipe.(nodename).pipeline.currgpuDevice == 0 || prmflow.pipe.(nodename).pipeline.kernellevel == 0
            % the node works on GPU#nextgpuDevice
            status.pipeline.(nodename).GPUdevice = prmflow.pipe.(nodename).pipeline.nextgpuDevice;
        else
            % the node works on GPU#currgpuDevice
            status.pipeline.(nodename).GPUdevice = prmflow.pipe.(nodename).pipeline.currgpuDevice;
        end
    else
        status.pipeline.(nodename).GPUdevice = 0;
    end

    % trace
    if status.debug.pooltrace_onoff && status.pipeline.(nodename).pipeline_onoff
        for jj = 1:length(dataflow.pipepool.(nodename))
            t = length(dataflow.pipepool.(nodename)(jj).trace);
            dataflow.pipepool.(nodename)(jj).trace(t+1).operator = 'prepare';
            dataflow.pipepool.(nodename)(jj).trace(t+1) = poolmirror(dataflow.pipepool.(nodename)(jj), ...
                dataflow.pipepool.(nodename)(jj).trace(t+1));
        end
    end
end

% call back the nodename
status.nodename = orignodename;

% to return the error message
if preperrorcode ~= 0
    status.errorcode = preperrorcode;
    status.errormsg = preperrormsg;
end

% done
status.jobdone = true;

end


function [carrynode, carryindex] = carryrecurse(pipepool, carrynode, carryindex)
% to find out the final carrypool

if pipepool.(carrynode)(carryindex).iscarried
    [carrynode, carryindex] = carryrecurse(pipepool, ...
        pipepool.(carrynode)(carryindex).carrynode, pipepool.(carrynode)(carryindex).carryindex);
end

end