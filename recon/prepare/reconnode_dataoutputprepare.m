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
    outputfiles = regexp(nodeprm.files, '(, +)|(,)', 'split');    
else
    % default output file is dicom image
    outputfiles = {'dicomimage'};
end
outputfiles_split = regexp(outputfiles, '_', 'split');
% name key
if ~isfield(nodeprm, 'namekey')
    if isfield(prmflow.protocol, 'namekey')
        namekey = ['_' prmflow.protocol.namekey];
    else
        namekey = '';
    end
end
% ext of corr files
if ~isfield(nodeprm, 'corrext')
    if isfield(prmflow.system, 'corrext')
        corrext = prmflow.system.corrext;
    else
        % default .corr
        corrext = '.corr';
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
% set original path
if strcmpi(prmflow.pipe.(nodename).outputpath, 'original')
    if ~isempty(prmflow.rawdata)
        prmflow.pipe.(nodename).outputpath = [fileparts(prmflow.rawdata) '\dicom\'];
    end
end
% mkdir outputpath
outputpath = prmflow.pipe.(nodename).outputpath;
if ~isfolder(outputpath)
    % mkdir if outputpath not exist
    mkdir(outputpath);
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

% default namerule
if isfield(nodeprm, 'namerule')
    namerule = nodeprm.namerule;
else
    namerule = 'manu';
end
% filename tags rule
filetagsrule = prmflow.system.filetagsrule;
% names set of the calibration tables
calinames = {'log2', 'air', 'beamharden', 'boneharden', 'nonlinear', 'crosstalk', 'offfocal', 'badchannel', 'hounsfield', ...
    'idealwater', 'detector'};
calinames = unique(cat(2, calinames, fieldnames(filetagsrule)'));
calinames = setdiff(calinames, {'raw', 'rawdata', 'image'});  % do not put rawdata and image in calibration set
% add 'corr' to the calinames
prmflow.pipe.(nodename).calinames = cellfun(@(x)[x 'corr'], calinames, 'UniformOutput', false);

Noutput = length(outputfiles);
% ini 
prmflow.pipe.(nodename).Noutput = Noutput;
prmflow.pipe.(nodename).OutputObject = struct();
prmflow.pipe.(nodename).OutputObject(Noutput) = struct();

dataflow.buffer.(nodename) = struct();
for ifile = 1 : Noutput
    prmflow.pipe.(nodename).OutputObject(ifile).outputfiles_split = outputfiles_split{ifile};
    obj_ifile = outputfiles_split{ifile}{1};
    % nametags
    if isfield(filetagsrule, obj_ifile)
        nametags = nametagrule(namerule, prmflow.protocol, filetagsrule.(obj_ifile));
    else
        nametags = nametagrule(namerule, prmflow.protocol);
    end
    switch lower(obj_ifile)
        case 'dicomimage'
            % I know the dicom dictionary has been set in reconinitial,
            % which was dicomdict('set', prmflow.system.dicomdictionary);
            % initial the the dcminfo in buffer
            dataflow.buffer.(nodename).dcminfo0 = getDicominfo(prmflow, status);
            dataflow.buffer.(nodename).dcmindex = 0;
            % file name
            prmflow.pipe.(nodename).OutputObject(ifile).filename = ...
                fullfile(outputpath, [outputfiles{ifile} namekey nametags]);
            prmflow.pipe.(nodename).OutputObject(ifile).fileext = '.dcm';
            prmflow.pipe.(nodename).OutputObject(ifile).fileversion = '';
            prmflow.pipe.(nodename).OutputObject(ifile).formatcfg = [];
            % object name
            prmflow.pipe.(nodename).OutputObject(ifile).name = 'dicomimage';
        case calinames
            % file name
            [filename, fileversion] = filerename(outputfiles{ifile}, outputfiles_split{ifile}, outputpath, ...
                nametags, corrext);
            prmflow.pipe.(nodename).OutputObject(ifile).filename = filename;
            prmflow.pipe.(nodename).OutputObject(ifile).fileext = corrext;
            prmflow.pipe.(nodename).OutputObject(ifile).fileversion = fileversion;
            % format configure
            prmflow.pipe.(nodename).OutputObject(ifile).formatcfg = readcfgfile(cfgmatchrule(filename, prmflow.IOstandard));
            % object name
            prmflow.pipe.(nodename).OutputObject(ifile).name = [lower(obj_ifile) 'corr'];
        case 'dataflow'
            if length(outputfiles_split{ifile}) > 1
                member = outputfiles_split{ifile}{2};
            else
                member = '';
            end
            % to same the dataflow.(member)
            switch member
                case {'rawdata', 'image'}
                    [filename, fileversion] = filerename(outputfiles{ifile}, outputfiles_split{ifile}, outputpath, ...
                        nametags, '.bin', 2);
                    prmflow.pipe.(nodename).OutputObject(ifile).filename = filename;
                    prmflow.pipe.(nodename).OutputObject(ifile).fileext = '.bin';
                    prmflow.pipe.(nodename).OutputObject(ifile).fileversion = fileversion;
                    % format configure
                    prmflow.pipe.(nodename).OutputObject(ifile).formatcfg = ...
                        readcfgfile(cfgmatchrule(filename, prmflow.IOstandard));
                    prmflow.pipe.(nodename).OutputObject(ifile).name = ['dataflow_', member];
                otherwise
                    % to same the dataflow or dataflow.(member) to .mat
                    filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '.mat']);
                    prmflow.pipe.(nodename).OutputObject(ifile).filename = filename;
                    prmflow.pipe.(nodename).OutputObject(ifile).fileext = '.mat';
                    prmflow.pipe.(nodename).OutputObject(ifile).fileversion = '';
                    prmflow.pipe.(nodename).OutputObject(ifile).formatcfg = [];
                    prmflow.pipe.(nodename).OutputObject(ifile).name = 'dataflow';
            end
        otherwise
            % to save any variable to .mat
            filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '.mat']);
            prmflow.pipe.(nodename).OutputObject(ifile).filename = filename;
            prmflow.pipe.(nodename).OutputObject(ifile).fileext = '.mat';
            prmflow.pipe.(nodename).OutputObject(ifile).fileversion = '';
            prmflow.pipe.(nodename).OutputObject(ifile).formatcfg = [];
            prmflow.pipe.(nodename).OutputObject(ifile).name = lower(obj_ifile);
    end
    % to overwrite the file
    prmflow.pipe.(nodename).OutputObject(ifile).wora = 'w';
    % 'a' is not supported yet
end

% pipe line
if pipeline_onoff
    % pipeline console paramters
    % carried
    if prmflow.protocol.tocarrythepools     % default is true
        prmflow.pipe.(nodename).pipeline.iscarried = true;
        % default was false
    end
else
    % while pipeline off,
    % default GPU off
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 0;
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function [filename, fileversion] = filerename(outputfiles, outputfiles_split, outputpath, nametags, ext, headnum)

if nargin < 6
    headnum = 1;
end

if (length(outputfiles_split) > headnum) && ...
        (regexp(outputfiles_split{end}, 'v\d+\.', 'once') == 1)
    % Umm, so we will name the file in this way #&!@@
    outfile_rename = linkstringcell(outputfiles_split(1:end-1), '_');
    % remove the first '_'
    outfile_rename = outfile_rename(2:end);
    outfile_rename(1) = upper(outfile_rename(1));
    outfile_rename = [outfile_rename nametags '_' outputfiles_split{end}];
    % done, that it is
    filename = fullfile(outputpath, [outfile_rename ext]);
    fileversion = outputfiles_split{end};
else
    % default version is 'v1.0'
    outfile_rename = outputfiles;
    outfile_rename(1) = upper(outfile_rename(1));
    filename = fullfile(outputpath, [outfile_rename nametags '_v1.0' ext]);
    fileversion = 'v1.0';
end

end

