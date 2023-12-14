function [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status)
% corr node, air correction prepare
% [dataflow, prmflow, status] = reconnode_airprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% parameters to use in prmflow
if ~isfield(nodeprm, 'blockwindow') || isempty(nodeprm.blockwindow)
    prmflow.raw.air.blockwindow = 10;
else
    prmflow.raw.air.blockwindow = nodeprm.blockwindow;
end
if pipeline_onoff
    prmflow.raw.air.blockskip = [0 0];
else
    prmflow.raw.air.blockskip = [prmflow.raw.air.blockwindow prmflow.raw.air.blockwindow];
end

% calibration table
aircorr = prmflow.corrtable.(nodename);
Nfocal = prmflow.raw.Nfocal;
Nref = aircorr.refnumber;
% reference error cut (ref-error-cut)
if isfield(aircorr, 'referrcut')
    referrcut = reshape(aircorr.referrcut, Nref, Nfocal);
else
    referrcut = zeros(Nref, Nfocal, 'single');
end
prmflow.raw.air.referrcut = single(referrcut);
% refpixel, refnumber
prmflow.raw.air.refpixel = single(aircorr.refpixel);
prmflow.raw.air.refnumber = single(aircorr.refnumber);

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    dataflow.buffer.(nodename) = struct();
    dataflow.buffer.(nodename).outpool = struct();
    dataflow.buffer.(nodename).ReadPoint = 1;
    dataflow.buffer.(nodename).WritePoint = 1;
    dataflow.buffer.(nodename).centernumber = 0;
    dataflow.buffer.(nodename).AvailPoint = 0;
    dataflow.buffer.(nodename).AvailViewindex = 0;
    dataflow.buffer.(nodename).poolsize = Inf;
    
    dataflow.buffer.(nodename).reflast = cell(1, Nfocal);
    switch prmflow.raw.scan
        case 'static'
            refblock = false(Nref, 0);
        case {'helical', 'halfaxial'}
            refblock = false(Nref, prmflow.raw.air.blockwindow);
        case 'axial'
            refblock = false(Nref, prmflow.raw.Nviewprot);
            dataflow.buffer.(nodename).poolsize = prmflow.raw.Nviewprot;
            dataflow.buffer.(nodename).WritePoint = prmflow.raw.Nviewprot - prmflow.raw.air.blockwindow + 1;
        otherwise
            % ??
            1;
    end
    dataflow.buffer.(nodename).outpool.rawhead.refblock = refblock;
    dataflow.buffer.(nodename).refblockrenewpoint = 1;
    dataflow.buffer.(nodename).Nrenew = 0;

    dataflow.buffer.(nodename).Axialpool = struct();
    dataflow.buffer.(nodename).currentview = 0;
end


% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end