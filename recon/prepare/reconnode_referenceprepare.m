function [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status)
% corr node, air correction prepare
% [dataflow, prmflow, status] = reconnode_referenceprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

if isfield(nodeprm, 'minviewnumber')
    prmflow.raw.air.minviewnumber = min(nodeprm.minviewnumber, prmflow.raw.Nviewprot/2);
else
    prmflow.raw.air.minviewnumber = min(200, prmflow.raw.Nviewprot/2);
end

% refernce error cut rescale
referrcut0 = prmflow.raw.air.referrcut;
if isfield(nodeprm, 'cutscale')
    prmflow.raw.air.referrcut = referrcut0.*nodeprm.cutscale;
else
    % default cut scale is 1.2;
    prmflow.raw.air.referrcut = referrcut0.*1.2;
end
% shall we scale it by mA?

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
refpixel = prmflow.raw.air.refpixel;
refpixelindex = [(1:refpixel) + pixelskip(1); (Npixel-refpixel+1:Npixel) - pixelskip(2)];
prmflow.raw.air.refpixelindex = refpixelindex;

% pipe line
if pipeline_onoff
    % initial input pool
    if ~isfield(dataflow.pipepool, nodename)
        dataflow.pipepool.(nodename) = status.defaultpooldata;
    end
    % initial private buffer
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).outpool = struct();
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = 1;
    dataflow.buffer.(nodename).AvailPoint = 0;
    dataflow.buffer.(nodename).AvailViewindex = 0;
    dataflow.buffer.(nodename).poolsize = Inf;
    
    dataflow.buffer.(nodename).reflast = cell(1, prmflow.raw.Nfocal);
    switch prmflow.raw.scan
        case 'static'
            refblock = false(prmflow.raw.air.refnumber, 0);
        case {'helical', 'halfaxial'}
            refblock = false(prmflow.raw.air.refnumber, prmflow.raw.air.blockwindow);
        case 'axial'
            refblock = false(prmflow.raw.air.refnumber, prmflow.raw.Nviewprot);
            dataflow.buffer.(nodename).poolsize = prmflow.raw.Nviewprot;
            dataflow.buffer.(nodename).WritePoint = prmflow.raw.Nviewprot - prmflow.raw.air.blockwindow + 1;
        otherwise
            % ??
            1;
    end
    dataflow.buffer.(nodename).outpool.rawhead.refblock = refblock;
    dataflow.buffer.(nodename).refblockrenewpoint = 1;
    dataflow.buffer.(nodename).Nrenew = 0;

    dataflow.buffer.(nodename).currentview = 0;
end


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end