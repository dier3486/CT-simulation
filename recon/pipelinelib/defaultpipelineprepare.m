function [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status)
% pipeline prepare lib function
% [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status);
% to fillup the prmflow.pipe.(nodename).pipeline by default settings and initial the dataflow.pipepool.(nodename).

% nodename
nodename = status.nodename;

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);

if ~status.pipeline.(nodename).pipeline_onoff
    % done while pipeline off
    return;
end

% to double (or int64)
prmflow.pipe.(nodename).pipeline = everything2single(prmflow.pipe.(nodename).pipeline, ...
    {'single', 'int32', 'int64', 'uint32', 'uint64'}, 'double');

% nodeprm
nodeprm = prmflow.pipe.(nodename);
nodepipeprm = nodeprm.pipeline;

% view rely strategy
if ischar(nodepipeprm.relystrategy) && strcmpi(nodepipeprm.relystrategy, 'NULL')
    prmflow.pipe.(nodename).pipeline.relystrategy = any(prmflow.pipe.(nodename).pipeline.viewrely);
    % The default relystrategy is 'None' (0) when none-viewrelying and 'Greedy'(1) when viewrelying;
end

% to 0 1 2
switch lower(prmflow.pipe.(nodename).pipeline.relystrategy)
    case {0, 'none'}
        prmflow.pipe.(nodename).pipeline.relystrategy = 0;
    case {1, 'greedy'}
        prmflow.pipe.(nodename).pipeline.relystrategy = 1;
    case {2, 'stingy'}
        prmflow.pipe.(nodename).pipeline.relystrategy = 2;
    otherwise
        error('Unknown view-rely relystrategy %s!', nodepipeprm.relystrategy);
end
% Normally, for axial the 'greedy viewrely' will ask the nextpool beging circulated and the 'stingy viewrely' will ask the
% currpool beging circulated; for helical the 'greedy viewrely' will write more data than the WritePoint be moved in nextpool
% and the 'stingy viewrely' will read more data than the ReadPoint be moved in currpool.

% viewcommon (the input views number must be an integer multiple of viewcommon)
if isnan(nodepipeprm.viewcommon)
    prmflow.pipe.(nodename).pipeline.viewcommon = lcm(prmflow.pipe.(nodename).pipeline.viewrescale(2), ...
        nodepipeprm.viewcommonfactor);
    % and multiple with Nfocal  
%     if isfield(prmflow, 'raw') && prmflow.raw.Nfocal > 0
%         prmflow.pipe.(nodename).pipeline.viewcommon = prmflow.pipe.(nodename).pipeline.viewcommon ...
%             * prmflow.raw.Nfocal;
%     end
    % don't do that, the Nfocal shall be employed in the node's prepare.
end
% The default viewcommon is the lcm(viewrescale(2), viewcommonfactor).
% Don't used to know those parameters while manually setting the viewcommon.

if mod(prmflow.pipe.(nodename).pipeline.viewrescale(2), prmflow.pipe.(nodename).pipeline.viewcommon) ~=0
    error('Unacceptalbe value of viewcommon. (You might to replace the viewcommon by viewcommonfactor.)');
end

% recover viewrely
viewrescale = prmflow.pipe.(nodename).pipeline.viewrescale;
if any(viewrescale~=1)
    prmflow.pipe.(nodename).pipeline.viewrely = ...
        ceil(prmflow.pipe.(nodename).pipeline.viewrely / viewrescale(2)) * viewrescale(2);
    prmflow.pipe.(nodename).pipeline.viewrely_out = ...
        prmflow.pipe.(nodename).pipeline.viewrely * viewrescale(1) / viewrescale(2);
else
    prmflow.pipe.(nodename).pipeline.viewrely_out = prmflow.pipe.(nodename).pipeline.viewrely;
end
% check viewextra/viewrely
if any(prmflow.pipe.(nodename).pipeline.viewextra > prmflow.pipe.(nodename).pipeline.viewrely_out)
    error('Too large viewextra! Which shall not more than the viewrely in output.');
end

% viewdelay is the views to be delayed by current node
prmflow.pipe.(nodename).pipeline.viewdelay = prmflow.pipe.(nodename).pipeline.viewdelay + ...
    prmflow.pipe.(nodename).pipeline.viewrely(2);

% output method
if strcmpi(nodepipeprm.outputmethod, 'NULL')
    if prmflow.pipe.(nodename).pipeline.relystrategy == 1
        % Greedy
        prmflow.pipe.(nodename).pipeline.outputmethod = 'cum';
    else
        prmflow.pipe.(nodename).pipeline.outputmethod = 'overwrite';
    end
end

% circulate mode by view-rely stragety
currcircshall = false;
nextcircshall = false;
if strcmpi(prmflow.protocol.scan, 'axial')
    if prmflow.pipe.(nodename).pipeline.viewdelay>0
        currcircshall = true;
    end
    if any(prmflow.pipe.(nodename).pipeline.viewrely > 0)
        switch prmflow.pipe.(nodename).pipeline.relystrategy
            case 1  % greedy
                % next pool should be circulate
                nextcircshall = true;
            case 2  % stingy
                % current pool should be circulate
                currcircshall = true;
            otherwise
                0;
        end
    end
end
% default currcirculte
if isempty(nodepipeprm.currcirculte)
    prmflow.pipe.(nodename).pipeline.currcirculte = currcircshall;
end
% default nextcirculte
if isempty(nodepipeprm.nextcirculte)
    prmflow.pipe.(nodename).pipeline.nextcirculte = nextcircshall;
end
% to check the asks by previous node
prevnode = status.pipeline.(nodename).prevnode;
Isprevexist = ~strcmpi(prevnode, 'NULL') && status.pipeline.(prevnode).pipeline_onoff;
if Isprevexist
    prevask = prmflow.pipe.(prevnode).pipeline;
else
    prevask = struct();
end
% circulte asking
if Isprevexist && ~isempty(prevask.nextcirculte)
    circulteask = prevask.nextcirculte;
else
    circulteask = false;
end
prmflow.pipe.(nodename).pipeline.currcirculte = prmflow.pipe.(nodename).pipeline.currcirculte | circulteask;
% viewdelayed
if Isprevexist
    prevviewrescale = prevask.viewrescale;
    prmflow.pipe.(nodename).pipeline.viewdelayed = ...
        floor((prevask.viewdelayed + prevask.viewdelay) * prevviewrescale(1) / prevviewrescale(2));
else
    prmflow.pipe.(nodename).pipeline.viewdelayed = 0;
end

% curr pool configured by user
dataflow.pipepool.(nodename)(1) = ...
    initialpool(dataflow.pipepool.(nodename)(1), prmflow.pipe.(nodename).pipeline, 'public');

% objecttype
if ~strcmpi(nodepipeprm.inputobjecttype, 'NULL')
    % user set
    objecttype = lower(nodepipeprm.inputobjecttype);
elseif Isprevexist && ~strcmpi(prevask.nextobjecttype, 'NULL')
    % objecttype asking
    objecttype = prevask.nextobjecttype;
    % I know if the user set the objecttype those asking will be ignored.
elseif Isprevexist && ~strcmpi(prevask.objecttype, 'NULL')
    % copy previous node's objecttype
    objecttype = prevask.objecttype;
else
    % default
    objecttype = 'rawdata';
end
prmflow.pipe.(nodename).pipeline.objecttype = objecttype;

% data fields
switch objecttype
    case 'rawdata'
        dataflow.pipepool.(nodename)(1).datafields = {'rawdata', 'rawhead'};
%         dataflow.pipepool.(nodename)(1).data = struct();
%         dataflow.pipepool.(nodename)(1).data.rawhead = struct();
%         dataflow.pipepool.(nodename)(1).data.rawdata = single([]);
    case 'image'
        dataflow.pipepool.(nodename)(1).datafields = {'image', 'imagehead'};
%         dataflow.pipepool.(nodename)(1).data = struct();
%         dataflow.pipepool.(nodename)(1).data.imagehead = struct();
%         dataflow.pipepool.(nodename)(1).data.image = single([]);
    case 'sparseimage'
        dataflow.pipepool.(nodename)(1).datafields = {'image', 'Sreconact', 'CSR', 'imagehead'};
%         dataflow.pipepool.(nodename)(1).data = struct();
%         dataflow.pipepool.(nodename)(1).data.imagehead = struct();
%         dataflow.pipepool.(nodename)(1).data.image = single([]);
%         dataflow.pipepool.(nodename)(1).data.Sreconact = logical([]);
    otherwise
        % user can set the datafields, but don't forget to set the objecttype to avoid those hard coded types.
        if isempty(dataflow.pipepool.(nodename)(1).datafields)
            dataflow.pipepool.(nodename)(1).datafields = {};
%             dataflow.pipepool.(nodename)(1).data = struct();
            % Note: empty datafields will forbid all the inputs.
        end
end

% to response the datasize asking for prm.pipeline
if Isprevexist && ~isnan(prevask.nextdatasize)
    prmflow.pipe.(nodename).pipeline.currdatasize = prevask.nextdatasize;
elseif isnan(prmflow.pipe.(nodename).pipeline.currdatasize)
    % we strongly suggest to define the currdatasize to claim the data in prepare
    if Isprevexist && ~isnan(prevask.currdatasize)
        prmflow.pipe.(nodename).pipeline.currdatasize = prevask.currdatasize;
    end
end

% default pipepool
if length(dataflow.pipepool.(nodename)) <= 1
    dataflow.pipepool.(nodename) = structmerge(dataflow.pipepool.(nodename), status.defaultpool, false, false);
end

% gpuDevice asking for prm.pipeline
if Isprevexist && prevask.nextgpuDevice > 0
    prmflow.pipe.(nodename).pipeline.currgpuDevice = prevask.nextgpuDevice;
end
% pool on GPU
if prmflow.pipe.(nodename).pipeline.currgpuDevice > 0
    dataflow.pipepool.(nodename)(1).bufferresource = sprintf('GPU device%d', prmflow.pipe.(nodename).pipeline.currgpuDevice);
    % I know the default value was 'CPU'
end
% Resource(type) asking for prm.pipeline
if Isprevexist && prevask.nextResource > 0
    prmflow.pipe.(nodename).pipeline.currResource = prevask.nextResource;
end
theResources = {'Matlab', 'sharedC'};
if prmflow.pipe.(nodename).pipeline.currResource > 0
    dataflow.pipepool.(nodename)(1).resourcetype = theResources{prmflow.pipe.(nodename).pipeline.currResource};
    % I know the default value was 'Matlab'
end

% currcirculte mode
if prmflow.pipe.(nodename).pipeline.currcirculte
    dataflow.pipepool.(nodename)(1).circulatemode = true;
end
% poolsize
if ~isavail(dataflow.pipepool.(nodename)(1).poolsize)
    if isavail(prmflow.pipe.(nodename).pipeline.currpoolsize)
        % user set
        dataflow.pipepool.(nodename)(1).poolsize = prmflow.pipe.(nodename).pipeline.currpoolsize;
    elseif Isprevexist && ~isnan(prevask.nextpoolsize)
        % poolsize asking for pipepool
        dataflow.pipepool.(nodename)(1).poolsize = prevask.nextpoolsize;
    elseif dataflow.pipepool.(nodename)(1).circulatemode
        % default circulated poolsize
        dataflow.pipepool.(nodename)(1).poolsize = prmflow.raw.viewpershot(1);
    else
        % default poolsize
        switch objecttype
            case 'rawdata'
                dataflow.pipepool.(nodename)(1).poolsize = prmflow.system.defaultrawpoolsize;
            case {'image', 'sparseimage'}
                dataflow.pipepool.(nodename)(1).poolsize = prmflow.system.defaultimagepoolsize;
            otherwise
                1;
        end
    end
end

% poolsize asking for pipepool
if Isprevexist && ~isnan(prevask.nextpoolsize)
    dataflow.pipepool.(nodename)(1).poolsize = max(dataflow.pipepool.(nodename)(1).poolsize, prevask.nextpoolsize);
end
% keepbottom asking for pipepool
if Isprevexist && ~isnan(prevask.nextkeepbottom)
    dataflow.pipepool.(nodename)(1).keepbottom = max(dataflow.pipepool.(nodename)(1).keepbottom, prevask.nextkeepbottom);
end
% basis asking for pipepool
if Isprevexist && prevask.nextdatabasis > 0
    dataflow.pipepool.(nodename)(1).databasis = prevask.nextdatabasis;
end

% basis asks to next node
if prmflow.pipe.(nodename).pipeline.nextdatabasis == 0
    % push on the databasis
    prmflow.pipe.(nodename).pipeline.nextdatabasis = dataflow.pipepool.(nodename)(1).databasis;
end

% device and resource asks to next node
if prmflow.pipe.(nodename).pipeline.GPUonoff || prmflow.pipe.(nodename).pipeline.CUDAonoff
    prmflow.pipe.(nodename).pipeline.nextgpuDevice = prmflow.system.GPUdevices(1);
    if prmflow.pipe.(nodename).pipeline.CUDAonoff
        % CUDA on
        prmflow.pipe.(nodename).pipeline.nextResource = 2;
    else
        prmflow.pipe.(nodename).pipeline.nextResource = 1;
    end
end

% check the carry by gpuDevice
if prmflow.pipe.(nodename).pipeline.currgpuDevice ~= prmflow.pipe.(nodename).pipeline.nextgpuDevice || ...
        prmflow.pipe.(nodename).pipeline.currResource ~= prmflow.pipe.(nodename).pipeline.nextResource || ...
        prmflow.pipe.(nodename).pipeline.nextdatabasis ~= dataflow.pipepool.(nodename)(1).databasis
    % can not carry
    prmflow.pipe.(nodename).pipeline.iscarried = false;
end

% carry asks (from prmflow..pipeline)
if prmflow.pipe.(nodename).pipeline.iscarried
    % ask circulte
    prmflow.pipe.(nodename).pipeline.nextcirculte = prmflow.pipe.(nodename).pipeline.currcirculte;
    % ask objecttype
    prmflow.pipe.(nodename).pipeline.nextobjecttype = prmflow.pipe.(nodename).pipeline.objecttype;
    % ask datasize
    if ~isnan(prevask.nextdatasize)
        prmflow.pipe.(nodename).pipeline.nextdatasize = prevask.nextdatasize;
    end
end

% carry asks (from dataflow..pipepool)
if prmflow.pipe.(nodename).pipeline.iscarried
    % ask poolsize
    if isnan(prmflow.pipe.(nodename).pipeline.nextpoolsize)
        prmflow.pipe.(nodename).pipeline.nextpoolsize = dataflow.pipepool.(nodename)(1).poolsize;
    end
    % ask keepbottom
    if dataflow.pipepool.(nodename)(1).keepbottom > 0 && isnan(prmflow.pipe.(nodename).pipeline.nextkeepbottom)
        prmflow.pipe.(nodename).pipeline.nextkeepbottom = dataflow.pipepool.(nodename)(1).keepbottom;
    end
    % not the GPU device and data basis can not forced setting by carry asking
%     % ask device (GPU)
%     if currGPUonoff
%         if prmflow.pipe.(nodename).pipeline.currgpuDevice > 1
%             prmflow.pipe.(nodename).pipeline.nextgpuDevice = prmflow.pipe.(nodename).pipeline.currgpuDevice;
%         else
%             prmflow.pipe.(nodename).pipeline.nextgpuDevice = prmflow.system.GPUdevices(1);
%         end
%     end
%     % ask basis
%     prmflow.pipe.(nodename).pipeline.nextdatabasis = dataflow.pipepool.(nodename)(1).databasis;
%     % but we suggest the prepare nodes should set those askings in not trivial cases.
end

% to initial the data in pipepool
if ~isnan(prmflow.pipe.(nodename).pipeline.currdatasize) && ~prmflow.pipe.(nodename).pipeline.iscarried
    datasize = prmflow.pipe.(nodename).pipeline.currdatasize;
    poolsize = dataflow.pipepool.(nodename)(1).poolsize;
    databasis = dataflow.pipepool.(nodename)(1).databasis;
    if poolsize>=0
        dataflow.pipepool.(nodename)(1).data = initialpooldata(dataflow.pipepool.(nodename)(1).data, objecttype, poolsize, datasize, databasis, ...
            prmflow.pipe.(nodename).pipeline.currgpuDevice, dataflow.pipepool.(nodename)(1).resourcetype);
    end
end

% debug
1;

end