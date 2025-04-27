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

flag_prepared = false;
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
    % to call the node's prepare function
    if any(exist(preparename) == [2 5 6]) % can run
        if echo_onoff
            if flag_prepared
                fprintf(', ');
            else
                fprintf(' [');
                flag_prepared = true;
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
    if status.pipeline.(status.nodename).pipeline_onoff
        [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status);
    end
end
if echo_onoff && flag_prepared
    fprintf(']');
end

% auto set carry nodes, only for pipepool(1) no branches
nodename = pipenodes{length(pipenodes)};
while ~strcmpi(nodename, 'NULL') && isfield(status.pipeline, nodename)
    % prevnode
    prevnode = status.pipeline.(nodename).prevnode;
    if strcmpi(prevnode, 'NULL')
        break;
    elseif ~status.pipeline.(prevnode).pipeline_onoff
        nodename = prevnode;
        continue;
    elseif ~isfield(prmflow.pipe, prevnode) || ~isfield(prmflow.pipe.(prevnode), 'pipeline')
        nodename = prevnode;
        continue;
    end
    if ~isfield(prmflow.pipe, nodename) || ~isfield(prmflow.pipe.(nodename), 'pipeline')
        nodename = prevnode;
        continue;
    end
    % is carry
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
    if ~status.pipeline.(nodename).pipeline_onoff
        continue;
    end
    for jj = 1:length(dataflow.pipepool.(nodename))
        if dataflow.pipepool.(nodename)(jj).iscarried
            carrynode = dataflow.pipepool.(nodename)(jj).carrynode;
            carryindex = dataflow.pipepool.(nodename)(jj).carryindex;
            if ~any(strcmp(dataflow.pipepool.(carrynode)(carryindex).carriages, nodename))
                dataflow.pipepool.(carrynode)(carryindex).carriages = ...
                    [dataflow.pipepool.(carrynode)(carryindex).carriages nodename];
                1;
            end
        end
    end
end

% double check the prepared nodes
for ii = 1:length(pipenodes)
    nodename = pipenodes{ii};
    if ~status.pipeline.(nodename).pipeline_onoff
        continue;
    end
    if isfield(prmflow.pipe, nodename) && isfield(prmflow.pipe.(nodename), 'pipeline')
        % check an illegal configure of the kernellevel 0
        if prmflow.pipe.(nodename).pipeline.kernellevel==0 && any(prmflow.pipe.(nodename).pipeline.viewrescale ~= 1) ...
                && prmflow.pipe.(nodename).pipeline.viewrescale(1) ~= prmflow.pipe.(nodename).pipeline.viewrescale(2)
            prmflow.pipe.(nodename).pipeline.kernellevel = 1;
            preperrorcode = 230;
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
            preperrorcode = 231;
            errormsg = ...
                sprintf('Node %s: illegal nodetype H-A.1.S.' ...
                , nodename);
            preperrormsg = [preperrormsg errormsg ' '];
        end
        % check illegal view extra
        if ~strcmpi(nextnode, 'NULL')
            if ~strcmpi(prmflow.pipe.(nodename).pipeline.nodetype, 'NULL') && status.pipeline.(nextnode).pipeline_onoff
                if any(prmflow.pipe.(nodename).pipeline.viewextra ~= 0) && (dataflow.pipepool.(nodename)(1).circulatemode || ...
                        dataflow.pipepool.(nextnode)(1).circulatemode)
                    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
                    preperrorcode = 232;
                    errormsg = ...
                        sprintf('Node %s, can not employ the view-extra in circulate-mode! Which is forcedly set to [0 0].' ...
                        , nodename);
                    preperrormsg = [preperrormsg errormsg ' '];
                end
            end
        end
        % check end node
        if strcmpi(nextnode, 'NULL') && prmflow.pipe.(nodename).pipeline.kernellevel==0 && ...
                status.pipeline.(nodename).pipeline_onoff
            % is it an error ?
%             preperrorcode = 239;
%             preperrormsg = sprintf('Node %s, TBC.' , pipenodes{ii});
        end
        % check carry node
%         if prmflow.pipe.(nodename).pipeline.iscarried && prmflow.pipe.(nodename).pipeline.kernellevel>0
%             preperrorcode = 233;
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
    end
end

% trace
if status.debug.pooltrace_onoff
    for ii = 1:length(pipenodes)
        nodename = pipenodes{ii};
        if ~status.pipeline.(nodename).pipeline_onoff
            continue;
        end
        if ~isfield(dataflow.pipepool, nodename)
            continue;
        end
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