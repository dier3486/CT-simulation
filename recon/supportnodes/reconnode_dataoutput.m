function [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status)
% support node, output data to file
% [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status);

% parameters set in pipe
outputprm = prmflow.pipe.(status.nodename);

% return if nothing to output
if ~isfield(outputprm, 'files') && isempty(outputprm.files)
    return;
end
% default namerule
if isfield(outputprm, 'namerule')
    namerule = outputprm.namerule;
else
    namerule = 'manu';
end
% filename tags rule
if isfield(prmflow.system, 'filetagsrule')
    filetagsrule = prmflow.system.filetagsrule;
    % I know the configure file of filetagsrule was read in reconinitial.m
else
    filetagsrule = struct();
end
% name key
if isfield(prmflow.protocol, 'namekey')
    namekey = ['_' prmflow.protocol.namekey];
else
    namekey = '';
end
% ext of corr files
if isfield(outputprm, 'corrext')
    corrext = outputprm.corrext;
elseif isfield(prmflow.system, 'corrext')
    corrext = prmflow.system.corrext;
else
    % default .corr
    corrext = '.corr';
end
% output path
if isfield(prmflow, 'outputpath')
    outputpath = prmflow.outputpath;
else
    % current path
    outputpath = '.';
end
if ~isfolder(outputpath)
    % mkdir if outputpath not exist
    mkdir(outputpath);
end
% IOstandard
if isfield(prmflow, 'IOstandard')
    IOstandard = prmflow.IOstandard;
else
    IOstandard = '';
end

% find out the outputs
outputfiles = regexp(outputprm.files, '(, +)|(,)', 'split');
outputfiles_split = regexp(outputfiles, '_', 'split');
% % name tags
% nametags = nametagrule(namerule, prmflow.protocol, [], prmflow.protocol.KV, prmflow.protocol.mA);

% ini the record of the output files
if ~isfield(prmflow, 'output')
    prmflow.output = struct();
end

% loop the outputs
for ifile = 1:length(outputfiles)
    outputobj = outputfiles_split{ifile}{1};
    % nametags
    if isfield(filetagsrule, outputobj)
        nametags = nametagrule(namerule, prmflow.protocol, filetagsrule.(outputobj));
    else
        nametags = nametagrule(namerule, prmflow.protocol);
    end
    switch lower(outputobj)
        case 'dicomimage'
            % to save diocom image(s)
            % I know the dicom dictionary has been set in reconinitial
%             if isfield(prmflow.system, 'dicomdictionary') && ~isempty(prmflow.system.dicomdictionary)
%                 dicomdict('set', prmflow.system.dicomdictionary);
%             end
            if isfield(dataflow, 'image')
                % dcminfo
                dcminfo = getDicominfo(prmflow, status);
                % image to output
                Nimage = size(dataflow.image, 3);
                prmflow.output.dicomimage = cell(Nimage, 1);
                for ii = 1:Nimage
                    % filename
                    filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags num2str(ii, '%03.0f') '.dcm']);
                    prmflow.output.dicomimage{ii} = filename;
                    % write dicom
                    image = uint16((dataflow.image(:,:,ii)'-1000-dcminfo(ii).RescaleIntercept)./dcminfo(ii).RescaleSlope);
                    dicomwrite(image, filename, dcminfo(ii), 'WritePrivate', true);
                end
            end
        case {'air', 'beamharden', 'boneharden', 'nonlinear', 'crosstalk', 'offfocal', 'badchannel', 'hounsefield', ...
                'idealwater', 'detector'}            
            % calibration tables
            objcorr = [outputobj 'corr'];
            if isfield(dataflow, objcorr)
                % filename
                if (length(outputfiles_split{ifile}) > 1) && (regexp(outputfiles_split{ifile}{end}, 'v\d+\.', 'once') ==1)
                    % Umm, so we will name the file in this way #&!@@
                    outfile_rename = linkstringcell(outputfiles_split{ifile}(1:end-1), '_');
                    outfile_rename = outfile_rename(2:end);
                    outfile_rename(1) = upper(outfile_rename(1));
                    outfile_rename = [outfile_rename nametags '_' outputfiles_split{ifile}{end}];
                    % done, that it is
                    filename = fullfile(outputpath, [outfile_rename corrext]);
                    fileversion = outputfiles_split{ifile}{end};
                else
                    % default version is 'v1.0'
                    outfile_rename = outputfiles{ifile};
                    outfile_rename(1) = upper(outfile_rename(1));
                    filename = fullfile(outputpath, [outfile_rename nametags '_v1.0' corrext]);
                    fileversion = 'v1.0';
                end
                % check version
                if ~checkfileversion(dataflow.(objcorr), fileversion)
                    warning('The corr ID of %s is not couple with the file format version %s!', objcorr, fileversion);
                end
                % pack the file
                prmflow.output.(objcorr) = filename;
                filecfg = readcfgfile(cfgmatchrule(filename, IOstandard));
                packstruct(dataflow.(objcorr), filecfg, filename);
            end
        case {'dataflow', 'prmflow', 'status'}
            % flow
            % save to .mat
            filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '.mat']);
            prmflow.output.(outputobj) = filename;
            save(filename, '-struct', outputobj);
        case {'raw', 'rawdata'}
            % pack the dataflow.rawdata and dataflow.rawhead back to .raw
            % not yet
            % TBC
        otherwise
            % save any variable in dataflow to .mat
            if isfield(dataflow, outputobj)
                filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '.mat']);
                prmflow.output.(outputobj) = filename;
                datatosave = dataflow.(outputobj);
                save(filename, 'datatosave');
            end
    end
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function r = checkfileversion(outputcorr, fileversion)
% check if the outputcorr.ID couple with fileversion

if isfield(outputcorr, 'ID')
    ID = outputcorr.ID(:);
else
    % no ID?
    r = true;
    return;
end

v = regexp(fileversion, '(\d+)(\.|$)', 'tokens');
v = cat(1, v{:});
v = cellfun(@str2num, v(:,1));

Nv = size(v, 1);
if Nv > 4
    r = false;
elseif any(ID(end-Nv+1:end)~=v)
    r = false;
else
    r = true;
end

end
