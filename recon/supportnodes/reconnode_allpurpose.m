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
nodeprm = prmflow.pipe.(nodename);

echo_onoff = nodeprm.pipeinfoecho;

if isfield(nodeprm, 'pipeline')
    pipeprm = nodeprm.pipeline;
else
    pipeprm = struct();
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% % is axial
% status.currentjob.isaxial = strcmpi(prmflow.raw.scan, 'axial');

if pipeline_onoff
    % next node
    nextnode = status.pipeline.(nodename).nextnode;
    % echo
    if echo_onoff
        if ~isempty(nodeprm.funct)
            fprintf(' (%s)', nodeprm.funct);
        end
        fprintf('\n        ');
    end

end

% copy pipeline_onoff to status.currentjob
% status.currentjob.pipeline_onoff = pipeline_onoff;

% prio
if pipeline_onoff
    if echo_onoff
        fprintf('PRIO: ');
    end

    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    
    if echo_onoff
        if pipeprm.priostep
            fprintf('read from previous %d views, ', status.currentjob.pipeline.readnumber);
        else
            fprintf('PRIO: skipped, ');
        end
    end
    
    if status.jobdone==0 || status.jobdone>=3
        % error or pass
        return;
    end
end

% node main
switch nodeprm.funct
    case {'loadrawdata', 'readraw'}
        [dataflow, prmflow, status] = fakenodereadrawdata(dataflow, prmflow, status);
    case 'donothing'
        status.jobdone = true;
        % I know the donothing node is level 0
    otherwise
        % default node function
        if pipeline_onoff
            % call kernel function
            switch nodeprm.pipeline.kernellevel
                case 0  % level 0 kernel function
                    0;
                    % nextnode -> nextnode
                    [dataflow.pipepool.(nextnode), prmflow, status] = ...
                        allpurposefakekernel(dataflow.pipepool.(nextnode), [], prmflow, status);

                case 1  % level 1 kernel function
                    1;
                    % currnode -> nextnode
                    [dataflow.pipepool.(nextnode), prmflow, status] = ...
                        allpurposefakekernel(dataflow.pipepool.(nextnode), dataflow.pipepool.(nodename), prmflow, status);
                otherwise
                    % error
            end
        else  % non-pipeline kernel function
            % dataflow->dataflow
            [dataflow, prmflow, status] = allpurposefakekernel([], dataflow, prmflow, status);
            status.jobdone = true;
        end
end

% post
if pipeline_onoff && pipeprm.poststep
    if echo_onoff
        % print the written view number
        fprintf('POST: write to next %d views and new available %d views, ', ...
            status.currentjob.pipeline.writenumber, status.currentjob.pipeline.newAvail);
    end

    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);

    if echo_onoff
        % print if write stuck
        if  ~isempty(dataflow.pipepool.(nextnode)) && dataflow.pipepool.(nextnode).WriteStuck
            fprintf('write locked ');
        end
    end

end

end