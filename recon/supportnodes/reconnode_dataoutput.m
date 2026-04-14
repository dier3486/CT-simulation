function [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status)
% support node, output data to file
% [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status);

% Copyright Dier Zhang
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

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
        calinames = nodeprm.calinames;
        Noutput = nodeprm.Noutput;
        OutputObject = nodeprm.OutputObject;
        % loop the outputs
        for ifile = 1 : Noutput
            % filename and format
            filename = OutputObject(ifile).filename;
            formatcfg = OutputObject(ifile).formatcfg;
            % output object
            outputobj = OutputObject(ifile).name;
            switch outputobj
                case 'dicomimage'
                    % to save diocom image(s)
                    % dcminfo
                    dcminfo0 = dataflow.buffer.(nodename).dcminfo0;
                    % imagesize
                    imagesize = prmflow.image.imagesize;
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
                        image_output = gather(dataflow.pipepool.(currnode)(currindex).data.image(:, index_output));
                    else
                        if ~isfield(dataflow, 'image')
                            continue;
                        elseif ~isfield(dataflow, 'imagehead')
                            warning('Miss the imagehead in dataflow while outputing the dicom images!');
                            dataflow.imagehead = struct();
                        end
                        Noutput = prmflow.image.Nimage;
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
                        filename_ii = [filename num2str(dcmindex + ii, '_%03.0f') '.dcm'];
                         % write dicom
                        dicomwrite(dcmimage, filename_ii, dcminfo(ii), 'WritePrivate', true);
                    end
                    dataflow.buffer.(nodename).dcmindex = dataflow.buffer.(nodename).dcmindex + Noutput;
                case calinames
                    % no calibration？
                    if ~isfield(dataflow, 'calibration')
                        continue;
                    end
                    % calibration tables
                    fileversion = OutputObject(ifile).fileversion;
                    % corr in dataflow
                    if isfield(dataflow.calibration, outputobj)
                        % check version
                        if ~checkfileversion(dataflow.calibration.(outputobj), fileversion)
                            warning('The corr ID of %s is not couple with the file format version %s!', outputobj, fileversion);
                        end
                        % pack the file
                        formatcfg = readcfgfile(cfgmatchrule(filename, prmflow.IOstandard));
                        packstruct(dataflow.calibration.(outputobj), formatcfg, filename);
                    end
                case {'prmflow', 'status'}
                    % save them to .mat
                    save(filename, '-struct', OutputObject(ifile).name);
                case {'dataflow'}
                    % save to .mat
                    if length(OutputObject(ifile).outputfiles_split) > 1
                        member = OutputObject(ifile).outputfiles_split{2};
                    else
                        member = '';
                    end
                    if isfield(dataflow, member)
                        save(filename, '-struct', dataflow.(member));
                    else
                        save(filename, '-struct', dataflow);
                    end
                case {'dataflow_rawdata', 'dataflow_image'}
                    if strcmp(OutputObject(ifile).name, 'dataflow_rawdata')
                        objectname = 'rawdata';
                        objecthead = 'rawhead';
                        headfields = prmflow.raw.rawheadfields;
                    elseif strcmp(OutputObject(ifile).name, 'dataflow_image')
                        objectname = 'image';
                        objecthead = 'imagehead';
                        headfields = [];
                        % shall be replaced by this
                        % headfields = prmflow.image.imageheadfields;
                    else
                        error('Wrong way of the case.');
                    end
                    % overwrite or add
                    wora = OutputObject(ifile).wora;
                    
                    % pack the dataflow.rawdata or pipepool.(node).data.rawdata and rawhead to .bin
                    if pipeline_onoff
                        plconsol = status.currentjob.pipeline;
                        % to output
                        IndexOutput = poolindex(dataflow.pipepool.(nodename)(1), plconsol.Index_in);
                        if dataflow.pipepool.(nodename)(1).iscarried
                            currnode = dataflow.pipepool.(nodename)(1).carrynode;
                            currindex = dataflow.pipepool.(nodename)(1).carryindex;
                        else
                            currnode = nodename;
                            currindex = 1;
                        end
                        % rawdata/image to write
                        datastruct = getdatastruct(dataflow.pipepool.(currnode)(currindex).data, IndexOutput, ...
                            objectname, objecthead, headfields, formatcfg);
                    else
                        if ~isfield(dataflow, 'image')
                            continue;
                        elseif ~isfield(dataflow, 'imagehead')
                            warning('Miss the imagehead in dataflow while outputing the dicom images!');
                            dataflow.imagehead = struct();
                        end
                        Noutput = prmflow.recon.Nimage;
                        IndexOutput = 1 : Noutput;
                        % rawdata/image to write
                        datastruct = getdatastruct(dataflow, IndexOutput, ...
                            objectname, objecthead, headfields, formatcfg);
                    end

                    % write to file
                    packstruct(datastruct, formatcfg, filename, wora);

                    if pipeline_onoff
                        % a+
                        prmflow.pipe.(nodename).OutputObject(ifile).wora = 'a';
                    end
                otherwise
                    0;
            end
        end

        % status
        if ~pipeline_onoff
            status.jobdone = true;
        end
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
elseif isfield(imagehead, 'ShotNumber')
    % the ShotNumber is AcquisitionNumber
    [dcminfo(:).AcquisitionNumber] = tac(imagehead.ShotNumber(:, index), 1);
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


function datastruct = getdatastruct(data, IndexOutput, objectname, objecthead, headfields, formatcfg)

if isempty(headfields)
    headfields = fieldnames(data.(objecthead))';
end
% headfields = setdiff(fieldnames(formatcfg), {'offset', 'class', 'size', 'number', objectname});

Nout = length(IndexOutput);
datastruct(Nout) = struct();

% rawdata/image basis
databasis = 1 + ~isreal(data.(objectname));     % 1 or 2
tagbasis = [objectname 'Basis'];
if isfield(formatcfg, tagbasis)
    [datastruct(:).(tagbasis)] = deal(databasis);
end
% There is not a class of complex, therefore the complex data can only be marked with the databasis.
% A complex single float will be formated in class='Single' and size=8.

% rawdata/image size
objectsize = classsize(lower(lower(formatcfg.(objectname).class)));
datasize = size(data.(objectname), 1) * objectsize;
% We shall re-check the datasize while databasis>2, now it will not happen.
tagsize = [objectname 'Size'];
if isfield(formatcfg, tagsize)       % shall always true
    [datastruct(:).(tagsize)] = deal(datasize);
else
    warning('Missing the data size in file format configure.');
end

headfields = setdiff(headfields, {tagbasis, tagsize});
% raw/image head
for ifield = headfields
    if isfield(formatcfg, ifield{1})
        fielddata = num2cell( gather(data.(objecthead).(ifield{1})(:, IndexOutput)), 1);
        [datastruct(:).(ifield{1})] = fielddata{:};
    end
end

% rawdata/image
fielddata = num2cell( gather(data.(objectname)(:, IndexOutput)), 1);
[datastruct(:).(objectname)] = fielddata{:};

end