function [dataflow, prmflow, status] = reconnode_datameanprepare(dataflow, prmflow, status)
% prepare node, prepare of data mean
% [dataflow, prmflow, status] = reconnode_datameanprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% nominalnumber (view number after mean)
if ~isfield(nodeprm, 'nominalnumber') || ~isavail(nodeprm.nominalnumber)
    prmflow.pipe.(nodename).nominalnumber = 1;
end

% skipviews
if ~isfield(nodeprm, 'viewskip') || isempty(nodeprm.viewskip)
    prmflow.pipe.(nodename).viewskip = 0;
end

% overshots (whether to do the data mean over the shots) 
if ~isfield(nodeprm, 'overshots') || isempty(nodeprm.overshots)
    prmflow.pipe.(nodename).overshots = false;
end

% recode the Nshot
prmflow.pipe.(nodename).Nshotbefore = prmflow.raw.Nshot;
prmflow.pipe.(nodename).viewpershotbefore = prmflow.raw.viewpershot;
prmflow.pipe.(nodename).startshotbefore = prmflow.raw.startshot;
prmflow.pipe.(nodename).Nviewbefore = prmflow.raw.Nview;

% after datamean
if prmflow.pipe.(nodename).overshots
    prmflow.raw.Nshot = 1;
    prmflow.raw.viewpershot = prmflow.pipe.(nodename).nominalnumber;
    prmflow.raw.startshot = 1;
    prmflow.raw.endshot = 1;
    prmflow.raw.Nview = prmflow.pipe.(nodename).nominalnumber;
    prmflow.raw.viewnumber = prmflow.pipe.(nodename).nominalnumber;
else
    prmflow.raw.viewpershot(prmflow.raw.startshot : prmflow.raw.endshot) = prmflow.pipe.(nodename).nominalnumber;
    prmflow.raw.Nview = prmflow.pipe.(nodename).nominalnumber * prmflow.raw.Nshot;
    prmflow.raw.viewnumber = prmflow.pipe.(nodename).nominalnumber * prmflow.raw.Nshot;

end
prmflow.raw.afterdatamean = true;
% Note: any nodes employed before the datamean shall pull those parameters to their private parameter set to use.

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % the datamean is A.1.N or H.1.N
    prmflow.pipe.(nodename).pipeline.kernellevel = 1;
    % prmflow.pipe.(nodename).pipeline.viewrescale = [0 1];
    prmflow.pipe.(nodename).pipeline.viewextra = [-prmflow.pipe.(nodename).viewskip, 0];
    % ask circulte 
    % if strcmpi(prmflow.raw.scan, 'axial')
        prmflow.pipe.(nodename).pipeline.nextcirculte = true;
    % end
    % ask poolsize 
    prmflow.pipe.(nodename).pipeline.nextpoolsize = prmflow.pipe.(nodename).nominalnumber;
    % no carry
    prmflow.pipe.(nodename).pipeline.iscarried = false;

    % default GPU off
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 0;
    end

    % dataflow.buffer.(nodename).Nstart = prmflow.pipe.(nodename).viewskip;

    dataflow.buffer.(nodename).ishot = 1;
    % dataflow.buffer.(nodename).Nsum = zeros(1, prmflow.pipe.(nodename).nominalnumber);
    
    
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end