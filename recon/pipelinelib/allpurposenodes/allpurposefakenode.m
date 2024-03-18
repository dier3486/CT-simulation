function [dataflow, prmflow, status] = allpurposefakenode(dataflow, prmflow, status, inputRinfo)

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% viewrely
if isfield(nodeprm, 'viewrely')
    viewrely = nodeprm.viewrely;
else
    viewrely = [0 0];
end

% end view
Nview = prmflow.raw.Nview;

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if pipeline_onoff
%     % did somthing read?
%     if inputRinfo.writenum == 0
%         % nothing to do
%         if inputRinfo.flagjobdone > 0
%             flagjobdone = 1;
%         else
%             % and nothing be ouput
%             flagjobdone = 3;
%         end
%         return;
%     end

    % next node
    nextnode = status.pipeline.(nodename).nextnode;
    if isempty(nextnode)
        nextnode = 'NULL';
    end

    % kernelInput
    kernelInput = struct();
    kernelInput.viewrely = viewrely;
    kernelInput.endviewindex = Nview;
    % but for images the endviewindex is not Nview. TBC

    status.jobdone = 0;
    nodelevel = nodeprm.level;
    switch nodelevel
        case 0
            % level 0
            kernelInput.AvailNumber = inputRinfo.writenum;
            kernelInput.StartPoint = inputRinfo.WritePoint;
            % the main data is in next pool
            [dataflow.pipepool.(nextnode), ~] = ...
                allpurposefakekernel(dataflow.pipepool.(nextnode), prmflow, status, kernelInput);
            % jobdone flag
            if kernelInput.AvailNumber > 0
                status.jobdone = 1;
            else
                % % nothing read
                status.jobdone = 3;
                % pass
            end
        case 1
            % level 1
            % new AvailNumber
            [kernelInput.StartPoint, kernelInput.AvailNumber] = ...
                newAvailNumber(dataflow.buffer.(nodename).outputpool, kernelInput.viewrely, kernelInput.endviewindex);
            % kernel function works on outputpool
            [dataflow.buffer.(nodename).outputpool.data, kernelRinfo] = ...
                allpurposefakekernel(dataflow.buffer.(nodename).outputpool.data, prmflow, status, kernelInput);
            % update outputpool.AvailNumber
            dataflow.buffer.(nodename).outputpool.AvailNumber = ...
                dataflow.buffer.(nodename).outputpool.AvailNumber + kernelRinfo.AvailNumber;
            % jobdone flag
            if kernelRinfo.AvailNumber > 0
                status.jobdone = 1; 
            else
                % nothing updated to output
                status.jobdone = 3;
                % pass
            end
        case 2
            % level 2
            % new AvailNumber
            [kernelInput.StartPoint, kernelInput.AvailNumber] = ...
                newAvailNumber(dataflow.buffer.(nodename).inputpool, kernelInput.viewrely, kernelInput.endviewindex);
            % kernel function works on inputpool
            [dataflow.buffer.(nodename).inputpool.data, kernelRinfo] = ...
                allpurposefakekernel(dataflow.buffer.(nodename).inputpool.data, prmflow, status, kernelInput);
            % update inputpool.AvailNumber
            dataflow.buffer.(nodename).inputpool.AvailNumber = ...
                dataflow.buffer.(nodename).inputpool.AvailNumber + kernelRinfo.AvailNumber;
            % copy data to output pool
            [dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).outputpool, ...
                dataflow.buffer.(nodename).outputpool.data, Rinfo] = ...
                pooldataoutput2(dataflow.buffer.(nodename).inputpool, dataflow.buffer.(nodename).inputpool.data, ...
                dataflow.buffer.(nodename).outputpool, dataflow.buffer.(nodename).outputpool.data);
            % update outputpool.AvailNumber
            dataflow.buffer.(nodename).outputpool.AvailNumber = ...
                dataflow.buffer.(nodename).outputpool.AvailNumber + Rinfo.writenum;
            % jobdone flag
            if Rinfo.writenum == 0
                % nothing updated to output
                status.jobdone = 3;
                % pass
            elseif Rinfo.writenum < Rinfo.datasize_towrite
                % some available data is stucked in inputpool
                status.jobdone = 2;
            else
                status.jobdone = 1;
            end
        otherwise

    end
else
    [dataflow, ~] = allpurposefakekernel(dataflow, prmflow, status);
    status.jobdone = true;
end


end


function [StartPoint, AvailNumber] = newAvailNumber(currpool, viewrely, endviewindex)
% new AvailNumber 

if nargin < 3
    endviewindex = inf;
end

if currpool.circulatemode
    % circulate mode
    RP2WP = poolp2p(currpool.ReadPoint, currpool.WritePoint, currpool.poolsize, 1);
    WE2WP = poolp2p(currpool.WriteEnd, currpool.WritePoint, currpool.poolsize, 1);
    if (WE2WP==1) && isfield(currpool, 'WriteStuck') && currpool.WriteStuck
        % axial data filled
        AvailNumber = RP2WP - currpool.AvailNumber + viewrely(1);
    else
        AvailNumber = RP2WP - currpool.AvailNumber - viewrely(2);
    end
    AvailNumber = max(AvailNumber, 0);
    StartPoint = mod(currpool.ReadPoint + currpool.AvailNumber - 1, currpool.poolsize) + 1;
else
    % not circulate
    if isfield(currpool, 'AvailViewindex')
        WV = currpool.AvailViewindex + currpool.WritePoint - currpool.AvailPoint;
    elseif isfield(currpool, 'ReadViewindex')
        WV = currpool.ReadViewindex + currpool.WritePoint - currpool.ReadPoint;
    else
        WV = 0;
    end
    if WV < endviewindex
        AP = currpool.WritePoint - viewrely(2) - 1;
    else
        AP = currpool.WritePoint - 1;
    end
    AvailNumber = AP - (currpool.ReadPoint + currpool.AvailNumber - 1);
    AvailNumber = max(AvailNumber, 0);
    StartPoint = currpool.ReadPoint + currpool.AvailNumber;
end

end


