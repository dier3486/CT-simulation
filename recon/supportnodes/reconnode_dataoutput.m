function [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status)
% support node, output data to file
% [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status);

% not prepared?
if ~status.pipeline.(status.nodename).prepared
    [dataflow, prmflow, status] = reconnode_dataoutputprepare(dataflow, prmflow, status);
    status.pipeline.(status.nodename).prepared = true;
end

% parameters set in pipe
nodename = status.nodename;

% pipeline_onoff
pipeline_onoff = status.currentjob.pipeline_onoff;

% prio
if pipeline_onoff
    % node prio-step
    [dataflow, prmflow, status] = nodepriostep(dataflow, prmflow, status);
    if status.currentjob.topass
        % error or pass
        return;
    end
end

% main
dataoutputKernelfuntion();

% post
if pipeline_onoff
    % post step
    [dataflow, prmflow, status] = nodepoststep(dataflow, prmflow, status);
end
% Done

% Kernel funtion
    function dataoutputKernelfuntion()
        % The anonymous function is static
        debug = [];

        % prepared parameters
        nodeprm = prmflow.pipe.(nodename);
        % outputfiles
        outputfiles = nodeprm.outputfiles;
        outputfiles_split = nodeprm.outputfiles_split;
        % namerule
        namerule = nodeprm.namerule;
        % name key
        namekey = nodeprm.namekey;

        % filename tags rule
        filetagsrule = prmflow.system.filetagsrule;

        % names set of the calibration tables
        calinames = {'air', 'beamharden', 'boneharden', 'nonlinear', 'crosstalk', 'offfocal', 'badchannel', 'hounsefield', ...
                        'idealwater', 'detector'};
        calinames = unique(cat(2, calinames, fieldnames(filetagsrule)'));
        calinames = setdiff(calinames, {'raw', 'rawdata'});
        
        % output path
        outputpath = nodeprm.outputpath;
        if ~isfolder(outputpath)
            % mkdir if outputpath not exist
            mkdir(outputpath);
        end
        % IOstandard
        IOstandard = prmflow.IOstandard;

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
                    % I know the dicom dictionary has been set in reconinitial,
                    % which was dicomdict('set', prmflow.system.dicomdictionary);
                    % dcminfo
                    if ~isfield(dataflow.buffer.(nodename), 'dcminfo')
                        dataflow.buffer.(nodename).dcminfo = getDicominfo(prmflow, status);
                    end
                    dcminfo0 = dataflow.buffer.(nodename).dcminfo;
                    % recon parameters
                    imagesize = prmflow.recon.imagesize;
                    
                    if pipeline_onoff
                        plconsol = status.currentjob.pipeline;
                        % to output
                        index_output = poolindex(dataflow.pipepool.(nodename)(1), plconsol.Index_in);
                        Noutput = length(index_output);
                        if dataflow.pipepool.(nodename)(1).iscarried
                            currnode = dataflow.pipepool.(nodename)(1).carrynode;
                            currindex = dataflow.pipepool.(nodename)(1).carryindex;
                        else
                            currnode = nodename;
                            currindex = 1;
                        end
                        % images to write
                        image_output = dataflow.pipepool.(currnode)(currindex).data.image(:, index_output);
                    else
                        if ~isfield(dataflow, 'image')
                            continue;
                        elseif ~isfield(dataflow, 'imagehead')
                            warning('Miss the imagehead in dataflow while outputing the dicom images!');
                            dataflow.imagehead = struct();
                        end
                        Noutput = prmflow.recon.Nimage;
                        index_output = 1 : Noutput;
                        % images to write
                        image_output = dataflow.image(:, index_output);
                    end
                    % with noise enhance
                    Isrealoutput = isreal(image_output);
                    
                    % dicom-info and images to write
                    if pipeline_onoff
                        dcminfo = getDicominfo2(dcminfo0, prmflow.recon, ...
                            dataflow.pipepool.(currnode)(currindex).data.imagehead, index_output, Isrealoutput);  
                    else
                        dcminfo = getDicominfo2(dcminfo0, prmflow.recon, ...
                            dataflow.imagehead, index_output, Isrealoutput);
                    end
                    % loop the output images
                    dcmindex = dataflow.buffer.(nodename).dcmindex;
                    for ii = 1 : Noutput
                        if Isrealoutput
                            dcmimage = uint16((image_output(:, ii) - 1000 - dcminfo0.RescaleIntercept) ...
                                ./ dcminfo0.RescaleSlope);
                        else
                            dcmimage = cat(2, ...
                                uint16((real(image_output(:, ii)) - 1000 - dcminfo0.RescaleIntercept) ./ dcminfo0.RescaleSlope), ...
                                uint16((imag(image_output(:, ii)) - dcminfo0.RescaleIntercept) ./ dcminfo0.RescaleSlope));
                        end
                        dcmimage = reshape(dcmimage, imagesize(2), imagesize(1), 1, []);
                        % filename
                        filename = fullfile(outputpath, ...
                            [outputfiles{ifile} namekey nametags num2str(dcmindex + ii, '_%03.0f') '.dcm']);
                         % write dicom
                        dicomwrite(dcmimage, filename, dcminfo(ii), 'WritePrivate', true);
                    end
                    dataflow.buffer.(nodename).dcmindex = dataflow.buffer.(nodename).dcmindex + Noutput;
                case calinames
                    % calibration tables
                    % ext of corr files
                    corrext = nodeprm.corrext;
                    % corr in dataflow
                    objcorr = [outputobj 'corr'];
                    if isfield(dataflow, objcorr)
                        % filename
                        if (length(outputfiles_split{ifile}) > 1) && ...
                                (regexp(outputfiles_split{ifile}{end}, 'v\d+\.', 'once') == 1)
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
                        if length(datatosave)==1 && isstruct(datatosave)
                            % to save the data in struct
                            save(filename, '-struct', 'datatosave');
                        else
                            save(filename, 'datatosave');
                        end
                    end
            end
        end

        % status
        status.jobdone = true;
    end

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


function dcminfo = getDicominfo2(dcminfo0, recon, imagehead, index, Isrealoutput)

if nargin < 5
    Isrealoutput = true;
end
if ~Isrealoutput
    % fix sopclassuid
    dcminfo0.sopclassuid = '1.2.840.10008.5.1.4.1.1.7.3';
end

Nimage = length(index);
dcminfo = repmat(dcminfo0, Nimage, 1);

% InstanceNumber
if isfield(imagehead, 'InstanceNumber')
    [dcminfo(:).InstanceNumber] = tac(imagehead.InstanceNumber(:, index), 1);
end
% AcquisitionNumber
if isfield(imagehead, 'AcquisitionNumber')
    [dcminfo(:).AcquisitionNumber] = tac(imagehead.AcquisitionNumber(:, index), 1);
elseif isfield(imagehead, 'Shot_Number')
    % the Shot_Number is AcquisitionNumber
    [dcminfo(:).AcquisitionNumber] = tac(imagehead.Shot_Number(:, index), 1);
end
% ImagePositionPatient
if isfield(imagehead, 'ImagePositionPatient')
    [dcminfo(:).ImagePositionPatient] = tac(imagehead.ImagePositionPatient(:, index), 1);
elseif isfield(imagehead, 'imagecenter')
    voxelsize = recon.voxelsize;
    XY = -voxelsize.*(recon.imagesize-1)./2;
    [dcminfo(:).ImagePositionPatient] = tac(imagehead.imagecenter(:, index) + [XY(:); 0], 1);
end
% SliceLocation
if isfield(imagehead, 'SliceLocation')
    [dcminfo(:).SliceLocation] = tac(imagehead.SliceLocation(:, index), 1);
elseif isfield(imagehead, 'imagecenter')
    [dcminfo(:).SliceLocation] = tac(imagehead.imagecenter(3, index) + recon.startcouch, 1);
end

end
