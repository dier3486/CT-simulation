function [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status)
% pipeline prepare lib function
% [dataflow, prmflow] = defaultpipelineprepare(dataflow, prmflow, status);
% to fillup the prmflow.pipe.(nodename).pipeline by default settings and initial the dataflow.pipepool.(nodename).

% nodename
nodename = status.nodename;

% fill empty
if ~isfield(prmflow.pipe.(nodename), 'pipeline') || isempty(prmflow.pipe.(nodename).pipeline)
    prmflow.pipe.(nodename).pipeline = struct();
end
% to double (or int64)
prmflow.pipe.(nodename).pipeline = everything2single(prmflow.pipe.(nodename).pipeline, ...
    {'single', 'int32', 'int64', 'uint32', 'uint64'}, 'double');

% nodeprm
nodeprm = prmflow.pipe.(nodename);
nodepipeprm = nodeprm.pipeline;
% these configure values will be put in prmflow.pipe.(nodename).pipeline
% prio/post step onoff
if ~isfield(nodepipeprm, 'priostep')
    if isempty(status.pipeline.(nodename).prevnode) || strcmpi(status.pipeline.(nodename).prevnode, 'NULL')
        prmflow.pipe.(nodename).pipeline.priostep = false;
        % defaultly skipping the prio-step for the first node.
    else
        prmflow.pipe.(nodename).pipeline.priostep = true;
    end
%     prmflow.pipe.(nodename).pipeline.priostep = true;
end
if ~isfield(nodepipeprm, 'poststep')
    prmflow.pipe.(nodename).pipeline.poststep = true;
end

% node level and kernel level
if isfield(nodepipeprm, 'kernellevel')
    kernellevel = nodepipeprm.kernellevel;
else
    kernellevel = 0;
end
prmflow.pipe.(nodename).pipeline.kernellevel = kernellevel;

% viewrely
if ~isfield(nodepipeprm, 'viewrely')
    prmflow.pipe.(nodename).pipeline.viewrely = [0 0];
end
% view rely strategy
if ~isfield(nodepipeprm, 'relystrategy')
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

% inputminlimit
if ~isfield(nodepipeprm, 'inputminlimit')
    prmflow.pipe.(nodename).pipeline.inputminlimit = 1;
end

% maxlimit
if ~isfield(nodepipeprm, 'inputmaxlimit')
    prmflow.pipe.(nodename).pipeline.inputmaxlimit = inf;
end

% viewdelay (only for axial)
if ~isfield(nodepipeprm, 'viewdelay')
    prmflow.pipe.(nodename).pipeline.viewdelay = 0;
end

% viewexpand (in output)
if ~isfield(nodepipeprm, 'viewexpand')
    prmflow.pipe.(nodename).pipeline.viewexpand = 0;
end

% viewextra (for helical head and bottom) which was rebin.Nextraview
if ~isfield(nodepipeprm, 'viewextra')
    prmflow.pipe.(nodename).pipeline.viewextra = [0 0];
end

% viewrescale (will change the viewnumber)
if ~isfield(nodepipeprm, 'viewrescale')
    prmflow.pipe.(nodename).pipeline.viewrescale = [1 1];
    % viewrescale = [Q(1) Q(2)], for that viewnum_out = viewnum_in*Q(1)/Q(2).
end

% viewcommon (the input views number must be an integer multiple of viewcommon)
if ~isfield(nodepipeprm, 'viewcommon') || isempty(nodepipeprm.viewcommon)
    if isfield(nodepipeprm, 'viewcommonfactor')
        prmflow.pipe.(nodename).pipeline.viewcommon = lcm(prmflow.pipe.(nodename).pipeline.viewrescale(2), ...
            nodepipeprm.viewcommonfactor);
    else
        prmflow.pipe.(nodename).pipeline.viewcommon = prmflow.pipe.(nodename).pipeline.viewrescale(2);
    end
    % and multiple with Nfocal
    if isfield(prmflow, 'raw')
        if prmflow.raw.Nfocal > 0
            prmflow.pipe.(nodename).pipeline.viewcommon = prmflow.pipe.(nodename).pipeline.viewcommon ...
                * prmflow.raw.Nfocal;
        end
    end
end
% The default viewcommon is the lcm(viewrescale(2), viewcommonfactor)*Nfocal.
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

% output method
if ~isfield(nodepipeprm, 'outputmethod')
    if prmflow.pipe.(nodename).pipeline.relystrategy == 1
        % Greedy
        prmflow.pipe.(nodename).pipeline.outputmethod = 'cum';
    else
        prmflow.pipe.(nodename).pipeline.outputmethod = 'overwrite';
    end
end

% % outputexpand (ask next node to maintain outputexpand more views' space in pool to get data)
% if ~isfield(nodepipeprm, 'outputexpand')
%     if prmflow.pipe.(nodename).pipeline.relystrategy == 1
%         prmflow.pipe.(nodename).pipeline.outputexpand = prmflow.pipe.(nodename).pipeline.viewrely(1);
%     else
%         prmflow.pipe.(nodename).pipeline.outputexpand = 0;
%     end
% end

% is a carry node
if ~isfield(nodepipeprm, 'iscarried')
    prmflow.pipe.(nodename).pipeline.iscarried = false;
    % will do this,
%     if defaultnode_onoff && kernellevel==0 && prmflow.pipe.(nodename).pipeline.relystrategy==0
%         prmflow.pipe.(nodename).pipeline.carrynode = true;
%     else
%         prmflow.pipe.(nodename).pipeline.carrynode = false;
%     end
end

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% initial pools
if pipeline_onoff
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
    if ~isfield(nodepipeprm, 'currcirculte')
        prmflow.pipe.(nodename).pipeline.currcirculte = currcircshall;
    end
    if ~isfield(nodepipeprm, 'nextcirculte')
        prmflow.pipe.(nodename).pipeline.nextcirculte = nextcircshall;
    end
    % to check the asks by previous node
    prevnode = status.pipeline.(nodename).prevnode;
    if strcmpi(prevnode, 'NULL')
        prevask = struct();
    elseif isfield(prmflow.pipe, prevnode) && isfield(prmflow.pipe.(prevnode), 'pipeline')
        prevask = prmflow.pipe.(prevnode).pipeline;
    else
        prevask = struct();
    end
    % circulte asking
    if isfield(prevask, 'nextcirculte')
        circulteask = prmflow.pipe.(prevnode).pipeline.nextcirculte;
    else
        circulteask = false;
    end
    prmflow.pipe.(nodename).pipeline.currcirculte = prmflow.pipe.(nodename).pipeline.currcirculte | circulteask;
    
    % initial the pipepool
    if ~isfield(dataflow.pipepool, nodename)
        dataflow.pipepool.(nodename) = struct();
    end
    % curr pool configured by user
    dataflow.pipepool.(nodename)(1) = ...
        initialpool(dataflow.pipepool.(nodename)(1), prmflow.pipe.(nodename).pipeline, 'public');
    
    % objecttype
    if isfield(nodepipeprm, 'inputobjecttype')
        % user set 
        objecttype = lower(nodepipeprm.inputobjecttype);
    elseif isfield(prevask, 'nextobjecttype')
        % objecttype asking
        objecttype = prevask.nextobjecttype;
        % I know if the user set the objecttype those asking will be ignored.
    elseif isfield(prevask, 'objecttype')
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
            dataflow.pipepool.(nodename)(1).data = struct();
            dataflow.pipepool.(nodename)(1).data.rawhead = struct();
            dataflow.pipepool.(nodename)(1).data.rawdata = single([]);   
        case 'image'
            dataflow.pipepool.(nodename)(1).datafields = {'image', 'imagehead'};
            dataflow.pipepool.(nodename)(1).data = struct();
            dataflow.pipepool.(nodename)(1).data.imagehead = struct();
            dataflow.pipepool.(nodename)(1).data.image = single([]);
        otherwise
            % user can set the datafields, but don't forget to set the objecttype to avoid those hard coded types.
            if ~isfield(dataflow.pipepool.(nodename)(1), 'datafields')
                dataflow.pipepool.(nodename)(1).datafields = {};
                dataflow.pipepool.(nodename)(1).data = struct();
            elseif isempty(dataflow.pipepool.(nodename)(1).datafields)
                dataflow.pipepool.(nodename)(1).datafields = {};
                dataflow.pipepool.(nodename)(1).data = struct();
                % Note: empty datafields will forbid all the inputs.
            end
    end

    % to response the datasize asking
    if isfield(prevask, 'nextdatasize')
        prmflow.pipe.(nodename).pipeline.currdatasize = prevask.nextdatasize;
    elseif ~isfield(prmflow.pipe.(nodename).pipeline, 'currdatasize')
        % we strongly suggest to define the currdatasize to claim the data in prepare
        if isfield(prevask, 'currdatasize')
            prmflow.pipe.(nodename).pipeline.currdatasize = prevask.currdatasize;
        end
    end

    % carry asks (from prmflow..pipeline)
    if prmflow.pipe.(nodename).pipeline.iscarried
        % ask circulte
        prmflow.pipe.(nodename).pipeline.nextcirculte = prmflow.pipe.(nodename).pipeline.currcirculte;
        % ask objecttype
        prmflow.pipe.(nodename).pipeline.nextobjecttype = prmflow.pipe.(nodename).pipeline.objecttype;
        % ask datasize
        if isfield(prevask, 'nextdatasize')
            prmflow.pipe.(nodename).pipeline.nextdatasize = prevask.nextdatasize;
        end
    end
    
    % default pipepool
    if length(dataflow.pipepool.(nodename)) <= 1
        dataflow.pipepool.(nodename) = structmerge(dataflow.pipepool.(nodename), status.defaultpool, false, false);
    end

    % poolsize and currcirculte mode
    if prmflow.pipe.(nodename).pipeline.currcirculte
        dataflow.pipepool.(nodename)(1).circulatemode = true;
        dataflow.pipepool.(nodename)(1).poolsize = prmflow.raw.viewpershot(1);
    elseif ~isavail(dataflow.pipepool.(nodename)(1).poolsize)
        % to set a default poolsize
        switch objecttype
            case 'rawdata'
                dataflow.pipepool.(nodename)(1).poolsize = prmflow.system.defaultrawpoolsize;
            case 'image'
                dataflow.pipepool.(nodename)(1).poolsize = prmflow.system.defaultimagepoolsize;
            otherwise
                1;
        end
    end

    % poolsize asking
    if isfield(prevask, 'nextpoolsize')
        dataflow.pipepool.(nodename)(1).poolsize = max(dataflow.pipepool.(nodename)(1).poolsize, prevask.nextpoolsize);
    end
    % keepbottom asking
    if isfield(prevask, 'nextkeepbottom')
        dataflow.pipepool.(nodename)(1).keepbottom = max(dataflow.pipepool.(nodename)(1).keepbottom, prevask.nextkeepbottom);
    end
    
    % carry asks (from dataflow.pipepool)
    if prmflow.pipe.(nodename).pipeline.iscarried
        % ask poolsize
        if ~isfield(prmflow.pipe.(nodename).pipeline, 'nextpoolsize')
            prmflow.pipe.(nodename).pipeline.nextpoolsize = dataflow.pipepool.(nodename)(1).poolsize;
        end
        % ask keepbottom
        if dataflow.pipepool.(nodename)(1).keepbottom > 0 && ~isfield(prmflow.pipe.(nodename).pipeline, 'nextkeepbottom')
            prmflow.pipe.(nodename).pipeline.nextkeepbottom = dataflow.pipepool.(nodename)(1).keepbottom;
        end
        % but we suggest the prepare nodes should set those askings in not trivial cases.
    end

    % to initial the data in pipepool
    if isfield(prmflow.pipe.(nodename).pipeline, 'currdatasize') && ~prmflow.pipe.(nodename).pipeline.iscarried
        datasize = prmflow.pipe.(nodename).pipeline.currdatasize;
        poolsize = dataflow.pipepool.(nodename)(1).poolsize;
        if datasize>=0 && poolsize>=0
            switch objecttype
                case 'rawdata'
                    dataflow.pipepool.(nodename)(1).data.rawdata = zeros(datasize, poolsize, 'single');
                case 'image'
                    dataflow.pipepool.(nodename)(1).data.image = zeros(datasize, poolsize, 'single');
                otherwise
                    1;
            end
        end
    end

end

