function [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status)
% corr node, air correction prepare
% [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if ~isfield(prmflow.correction, 'air')
    % did not do the air correction?
    status.jobdone = false;
    status.errorcode = 220;
    status.errormsg = 'The reference correction can not be employed before or without air correction.';
    return;
end

% refernce error cut rescale
referrcut0 = prmflow.correction.air.referrcut;
if isfield(nodeprm, 'cutscale')
    prmflow.correction.air.referrcut = referrcut0.*nodeprm.cutscale;
else
    % default cut scale is 1.2;
    prmflow.correction.air.referrcut = referrcut0.*1.2;
end
% shall we scale it by mA?

% over write blockwindow
if strcmpi(prmflow.raw.scan, 'static')
    % ignore the blockwindow for static scanning
    prmflow.correction.air.blockwindow = 0;
else
    if isfield(nodeprm, 'blockwindow')
        prmflow.correction.air.blockwindow = nodeprm.blockwindow;
    elseif ~isfield(prmflow.correction.air, 'blockwindow')
        prmflow.correction.air.blockwindow = 12;
    end
end

% refernceindex
if isfield(nodeprm, 'refpixelskip')
    pixelskip = nodeprm.refpixelskip;
    if length(pixelskip)==1
        pixelskip = [pixelskip pixelskip];
    end
else
    pixelskip = [0 0];
end
Npixel = prmflow.raw.Npixel;
refpixel = prmflow.correction.air.refpixel;
refpixelindex = [(1:refpixel) + pixelskip(1); (Npixel-refpixel+1:Npixel) - pixelskip(2)];
prmflow.correction.air.refpixelindex = refpixelindex;

% slice independent refernece?
if isfield(nodeprm, 'sliceindependent')
    prmflow.correction.air.sliceindependent = nodeprm.sliceindependent;
else
    prmflow.correction.air.sliceindependent = false;
end

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % the reference correction is H-H.0.G or A.0.G 
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    if strcmpi(prmflow.protocol.scan, 'static')
        % but static scan is in type is H.0.N specially
        prmflow.pipe.(nodename).pipeline.viewrely = [0 0];
        prmflow.pipe.(nodename).pipeline.relystrategy = 0;
    else
        prmflow.pipe.(nodename).pipeline.viewrely = [prmflow.correction.air.blockwindow prmflow.correction.air.blockwindow];
        prmflow.pipe.(nodename).pipeline.relystrategy = 'greedy';
        prmflow.pipe.(nodename).pipeline.outputmethod = 'cum';
    end
    if strcmpi(prmflow.protocol.scan, 'axial')
        prmflow.pipe.(nodename).pipeline.nextcirculte = true;
    end

    % inputminlimit (nodeprm.minviewnumber)
    if isfield(nodeprm, 'minviewnumber')
        prmflow.pipe.(nodename).pipeline.inputminlimit = min(nodeprm.minviewnumber, prmflow.raw.Nviewprot/2);
    else
        prmflow.pipe.(nodename).pipeline.inputminlimit = min(200, prmflow.raw.Nviewprot/2);
    end
    

    % private buffer (to save the reflast, only)
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).reflast = cell(1, prmflow.raw.Nfocal);
end

end