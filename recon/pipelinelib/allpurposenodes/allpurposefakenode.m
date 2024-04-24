function [dataflow, prmflow, status] = allpurposefakenode(dataflow, prmflow, status, inputRinfo)

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% kernelInput
kernelInput = struct();
% viewrely
if isfield(nodeprm, 'viewrely')
    kernelInput.viewrely = nodeprm.viewrely;
else
    kernelInput.viewrely = [0 0];
end
% view rely strategy
if isfield(nodeprm, 'relystrategy')
    kernelInput.relystrategy = nodeprm.relystrategy;
    % relystrategy = 'Greedy'(1) , 'Stingy'(2) or 'Null'(0)
else
    kernelInput.relystrategy = any(kernelInput.viewrely);
end
% extraview 
% not now, TBC
% object type
if isfield(nodeprm, 'objecttype')
    kernelInput.objecttype = nodeprm.objecttype;
    % objecttype = 'Rawdata'(1) or 'Image'(2).
else
    kernelInput.objecttype = 'Rawdata';
end
% endviewindex
switch lower(kernelInput.objecttype)
    case {'rawdata', 1}
        kernelInput.endviewindex = prmflow.raw.Nview;
    case {'image', 2}
        kernelInput.endviewindex = prmflow.recon.Nimage;
    otherwise
        kernelInput.endviewindex = inf;
end

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

    status.jobdone = 0;
    nodelevel = nodeprm.level;
    switch nodelevel
        case 0
            % level 0
            kernelInput.DataStart = inputRinfo.origWritePoint;
            kernelInput.DataEnd = kernelInput.DataStart + inputRinfo.writenum;
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
            kernelInput = ...
                newAvailNumber(kernelInput, dataflow.buffer.(nodename).outputpool);
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
            kernelInput = ...
                newAvailNumber(kernelInput, dataflow.buffer.(nodename).inputpool);
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


function kernelInput = newAvailNumber(kernelInput, currpool)
% new AvailNumber 

if currpool.circulatemode
    % circulate mode
    RP2WP = poolp2p(currpool.ReadPoint, currpool.WritePoint, currpool.poolsize, 1);
    WE2WP = poolp2p(currpool.WriteEnd, currpool.WritePoint, currpool.poolsize, 1);
    if (WE2WP==1) && isfield(currpool, 'WriteStuck') && currpool.WriteStuck
        % axial data filled
        AvailNumber = RP2WP - currpool.AvailNumber + kernelInput.viewrely(1);
    else
        AvailNumber = RP2WP - currpool.AvailNumber - kernelInput.viewrely(2);
    end
    kernelInput.AvailNumber = max(AvailNumber, 0);
    kernelInput.DataStart = mod(currpool.ReadPoint + currpool.AvailNumber - 1, currpool.poolsize) + 1;
    kernelInput.DataEnd = kernelInput.DataStart + kernelInput.AvailNumber;
    % rely strategy TBC 
else
    % not circulate
    if isfield(currpool, 'AvailViewindex')
        WV = currpool.AvailViewindex + currpool.WritePoint - currpool.AvailPoint;
    elseif isfield(currpool, 'ReadViewindex')
        WV = currpool.ReadViewindex + currpool.WritePoint - currpool.ReadPoint;
    else
        WV = 0;
    end
    if WV < kernelInput.endviewindex
        AP = currpool.WritePoint - viewrely(2) - 1;
    else
        AP = currpool.WritePoint - 1;
    end
    AvailNumber = AP - (currpool.ReadPoint + currpool.AvailNumber - 1);
    kernelInput.AvailNumber = max(AvailNumber, 0);
    kernelInput.DataStart = currpool.ReadPoint + currpool.AvailNumber;
    kernelInput.DataEnd = kernelInput.DataStart + kernelInput.AvailNumber;
    % rely strategy TBC 
end

end


function [dataflow, Rinfo] = allpurposefakekernel(dataflow, prmflow, status, Ininfo)
% a fake node kernel function

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% operator
if isfield(nodeprm, 'mainopt')
    mainopt = prmflow.pipe.(nodename).mainopt;
else
    mainopt = @(x) x;
end

if nargin<4
    if isfield(dataflow, 'rawdata')
        dataflow.rawdata = mainopt(dataflow.rawdata);
    end
    if isfield(dataflow, 'image')
        dataflow.image = mainopt(dataflow.image);
    end
    Rinfo = [];
else
    DataStart = Ininfo.DataStart;
    DataEnd = Ininfo.DataEnd;

    if isfield(dataflow, 'rawdata')
        dataflow.rawdata(:, DataStart:DataEnd) = mainopt(dataflow.rawdata(:, DataStart:DataEnd));
        % not support circulation
    end
%     if isfield(dataflow, 'image') && (nargin>2) && ~isempty(imagestartend)
%         dataflow.image(:, imagestartend(1):imagestartend(2)) = dataflow.image(:, imagestartend(1):imagestartend(2)) + 1;
%     end
    % return Rinfo, ummm seems nothing to return
    Rinfo = Ininfo;
end

end

