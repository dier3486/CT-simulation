function [dataflow, prmflow, status] = reconnode_allpurposeprepare(dataflow, prmflow, status)
% prepare node, prepare to do nothing
% [dataflow, prmflow, status] = reconnode_allpurposeprepare(dataflow, prmflow, status);
% wow~

% parameters set in pipe
nodename = status.nodename;
if isfield(prmflow.pipe, nodename)
    nodeprm = prmflow.pipe.(nodename);
else
    nodeprm = struct();
end

% main operator
if isfield(nodeprm, 'mainopt')
    if ischar(nodeprm.mainopt)
        prmflow.pipe.(nodename).mainopt = eval(nodeprm.mainopt);
    end
else
    prmflow.pipe.(nodename).mainopt = @(x) x+1;
end

% % objecttype
% if ~isfield(nodeprm, 'inputobjecttype')
%     prmflow.pipe.(nodename).inputobjecttype = 'rawdata';
% end

% echo onoff
if ~isfield(nodeprm, 'pipeinfoecho')
    prmflow.pipe.(nodename).pipeinfoecho = true;
end

% node function
if isfield(nodeprm, 'funct')
    nodefunct = nodeprm.funct;
else
    prmflow.pipe.(nodename).funct = '';
    nodefunct = '';
end

% prepare by existing nodes
if isempty(nodefunct)
    preparename = '';
else
    preparename = ['reconnode_' lower(nodefunct) 'prepare'];
end
if any(exist(preparename) == [2 5 6]) % can run
    if status.debug.echo_onoff
        fprintf(['(' nodefunct ')']);
    end
    preparefun = str2func(preparename);
    [dataflow, prmflow, status] = preparefun(dataflow, prmflow, status);
    if status.jobdone==0
        return;
    end
elseif isfield(nodeprm, 'pipeline') && isfield(prmflow, 'prepare')
    % objecttype
    if ~isfield(nodeprm.pipeline, 'inputobjecttype')
        prmflow.pipe.(nodename).pipeline.inputobjecttype = 'rawdata';
    end
    % a default prepare of prmflow.prepare
    switch lower(prmflow.pipe.(nodename).pipeline.inputobjecttype)
        case 'rawdata'
            if isfield(nodeprm.pipeline, 'viewrescale') && any(nodeprm.pipeline.viewrescale ~= 1)
                Q = nodeprm.pipeline.viewrescale;
                prmflow.prepare.viewpershot = prmflow.prepare.viewpershot.*Q(1)./Q(2);
                prmflow.prepare.Nview = prmflow.prepare.Nview.*Q(1)./Q(2);
            end
            if isfield(nodeprm.pipeline, 'viewextra') && any(nodeprm.pipeline.viewextra ~= 0)
                E = nodeprm.pipeline.viewextra;
                prmflow.prepare.viewpershot = prmflow.prepare.viewpershot + E(1) + E(2);
                prmflow.prepare.Nview = prmflow.prepare.Nview + (E(1) + E(2))*prmflow.prepare.Nshot;
            end
        case 'image'
            2;
            % TBC
        otherwise
            0;
            % TBC
    end
end


status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


