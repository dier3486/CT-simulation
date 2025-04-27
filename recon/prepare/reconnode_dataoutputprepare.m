function [dataflow, prmflow, status] = reconnode_dataoutputprepare(dataflow, prmflow, status)
% prepare node, dataoutput
% [dataflow, prmflow, status] = reconnode_dataoutputprepare(dataflow, prmflow, status);

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);
% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prepare of dataoutput
% outputfiles
if isfield(nodeprm, 'files')
    prmflow.pipe.(nodename).outputfiles = regexp(nodeprm.files, '(, +)|(,)', 'split');    
else
    % default output file is dicom image
    prmflow.pipe.(nodename).outputfiles = {'dicomimage'};
end
prmflow.pipe.(nodename).outputfiles_split = regexp(prmflow.pipe.(nodename).outputfiles, '_', 'split');
% default namerule
if ~isfield(nodeprm, 'namerule')
    prmflow.pipe.(nodename).namerule = 'manu';
end
% name key
if ~isfield(nodeprm, 'namekey')
    if isfield(prmflow.protocol, 'namekey')
        prmflow.pipe.(nodename).namekey = ['_' prmflow.protocol.namekey];
    else
        prmflow.pipe.(nodename).namekey = '';
    end
end
% ext of corr files
if ~isfield(nodeprm, 'corrext')
    if isfield(prmflow.system, 'corrext')
        prmflow.pipe.(nodename).corrext = prmflow.system.corrext;
    else
        % default .corr
        prmflow.pipe.(nodename).corrext = '.corr';
    end
end
% output path
if ~isfield(nodeprm, 'outputpath')
    if isfield(prmflow.protocol, 'outputpath')
        prmflow.pipe.(nodename).outputpath = prmflow.protocol.outputpath;
    else
        % current path
        prmflow.pipe.(nodename).outputpath = '.';
    end
end
% set orign path
if strcmpi(prmflow.pipe.(nodename).outputpath, 'orign')
    if ~isempty(prmflow.rawdata)
        outputpath = [fileparts(prmflow.rawdata) '\dicom\'];
        if ~isfolder(outputpath)
            mkdir(outputpath);
        end
        prmflow.pipe.(nodename).outputpath = outputpath;
    end
end

% outputblock (imageblock)
if ~isfield(nodeprm, 'outputblock') || isempty(nodeprm.outputblock)
    prmflow.pipe.(nodename).outputblock = 1;
end
outputblock = prmflow.pipe.(nodename).outputblock;

if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    prmflow.pipe.(nodename).pipeline.kernellevel = 0;
    % minoutput
    if ~isfield(prmflow.pipe.(nodename).pipeline, 'outputminlimit')
        prmflow.pipe.(nodename).pipeline.outputminlimit = outputblock;
    end
    % rescale
    prmflow.pipe.(nodename).pipeline.viewrescale = [outputblock outputblock];
end

dataflow.buffer.(nodename) = struct();
% dicom
if any(strcmpi([prmflow.pipe.(nodename).outputfiles_split{:}], 'dicomimage'))
    dataflow.buffer.(nodename).dcminfo = getDicominfo(prmflow, status);
    dataflow.buffer.(nodename).dcmindex = 0;
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end
