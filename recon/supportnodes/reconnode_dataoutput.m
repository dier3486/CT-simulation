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
    namerule = '';
end
% name key
if isfield(prmflow.protocol, 'namekey')
    namekey = ['_' prmflow.protocol.namekey];
else
    namekey = '';
end
% output path
if isfield(prmflow, 'outputpath')
    outputpath = prmflow.outputpath;
else
    % current path
    outputpath = '.';
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
% name tags
nametags = nametagrule(namerule, prmflow.protocol, [], prmflow.protocol.KV, prmflow.protocol.mA);

% ini the record of the output files
if ~isfield(prmflow, 'output')
    prmflow.output = struct();
end

% loop the outputs
for ifile = 1:length(outputfiles)
    outputobj = outputfiles_split{ifile}{1};
    switch lower(outputobj)
        case 'dicomimage'
            % diocom image
            if isfield(dataflow, 'dicomimage')
                if iscell(dataflow.image)
                    % to record the output files
                    prmflow.output.dicomimage = cell(1, length(dataflow.image));
                    % loop the images
                    for ii = 1:length(dataflow.image)
                        % image to output
                        image = dataflow.image{ii};
                        image = reshape(image, size(image, 1), []);
                        if ~isinteger(image)
                            image = int16(image);
                        end
                        % filename
                        filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '_' num2str(ii) '.dcm']);
                        prmflow.output.dicomimage{ii} = filename;
                        % diocom
                        dicomwrite(image, filename);
                        % dicominfo is not supported yet, TBC
                    end
                else
                    % image to output
                    image = dataflow.image;
                    image = reshape(image, size(image, 1), []);
                    if ~isinteger(image)
                        image = int16(image);
                    end
                    % filename
                    filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '.dcm']);
                    prmflow.output.dicomimage = filename;
                    % diocom
                    dicomwrite(image, filename);
                    % dicominfo is not supported yet, TBC
                end  
            end
        case {'air', 'beamharden', 'boneharden', 'nonlinear', 'crosstalk', 'offfocal', 'badchannel', 'housefield'}
            % calibration tables
            objcorr = [outputobj 'corr'];
            if isfield(dataflow, objcorr)
                % filename
                if (length(outputfiles_split{ifile}) > 1) && (regexp(outputfiles_split{ifile}{end}, 'v\d+\.', 'once') ==1)
                    % Umm, so we will name the file in this way #&!@@
                    Nsplit = length(outputfiles_split{ifile});
                    outfile_rename = cell(1, Nsplit*2-3);
                    outfile_rename(1:2:(Nsplit*2-3)) = outputfiles_split{ifile}(1:end-1);
                    if Nsplit>2
                        outfile_rename{2:2:(Nsplit*2-4)} = {'_'};
                    end
                    outfile_rename = [outfile_rename{:} namekey nametags '_' outputfiles_split{ifile}{end}];
                    % done, that it is
                    filename = fullfile(outputpath, [outfile_rename '.corr']);
                else
                    % default version is 'v1.0'
                    filename = fullfile(outputpath, [outputfiles{ifile} namekey nametags '_v1.0' '.corr']);
                end
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
        case 'raw'
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

