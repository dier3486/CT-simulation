function [dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status)
% prepare node, BP prepare, set FOV, image size, center (XYZ), tilt ... for the images
% [dataflow, prmflow, status] = reconnode_backprojectionprepare(dataflow, prmflow, status);

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

% parameters set in pipe
nodename = status.nodename;
nodeprm = prmflow.pipe.(nodename);

% pipeline_onoff
pipeline_onoff = status.pipeline.(nodename).pipeline_onoff;

% status ini
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

% scan
scan = lower(prmflow.protocol.scan);
% recon method
if isfield(nodeprm, 'method') && ~isempty(nodeprm.method)
    recon_method = nodeprm.method;
else
    % default BP method
    switch scan
        case 'axial'
            recon_method = '2D';
        case 'helical'
            recon_method = '';
        case 'conveyor'
            recon_method = '';
        otherwise
            recon_method = '';
    end
end

if ~strncmpi(recon_method, scan, length(scan))
    prmflow.recon.method = [scan recon_method];
else
    prmflow.recon.method = recon_method;
end
% The default recon method of axial is Axial2D

% common prepare
prmflow.recon = commonbpprepare(prmflow.recon, prmflow.protocol, prmflow.system, nodeprm);

% switch recon method
switch lower(prmflow.recon.method)
    case {'axial2d', 'axial3d'}
        prmflow.recon.methodswitch = 1;
        prmflow.recon = axial3Dprepare(prmflow.recon, nodeprm);
    case {'axialhalf'}
        prmflow.recon.methodswitch = 2;
        % TBC
    case {'helical', 'helical3d', 'helicalpiline', 'helicalsuper', 'helicalpilinesuper' , ...
            'helicalbigpitch'}
        prmflow.recon.methodswitch = 3;
        % helical is always 3D
        prmflow.recon = helicalprepare(prmflow.recon, nodeprm);
    
    case {'conveyor', 'conveyor3d', 'conveyorpiline'}
        prmflow.recon.methodswitch = 4;
        % conveyor(piline)
        prmflow.recon = conveyorprepare(prmflow.recon, nodeprm);
        % ask TableGear
        if ~isincell(prmflow.raw.rawheadfields, 'TableEncoder')
            prmflow.raw.rawheadfields = [prmflow.raw.rawheadfields, {'TableEncoder'}];
        end
        if ~isincell(prmflow.raw.rawheadfields, 'TableGear')
            prmflow.raw.rawheadfields = [prmflow.raw.rawheadfields, {'TableGear'}];
        end
    otherwise
        prmflow.recon.methodswitch = 0;
        % do nothing
        1;
        % no topo
end

if contains(lower(prmflow.recon.method), 'super')
    prmflow.pipe.(nodename).operator = 'super';
    if prmflow.recon.previewed
        prmflow.recon.Spresuper = ~prmflow.recon.Spreview & prmflow.recon.Ssuper(:, 1);
    end
else
    prmflow.pipe.(nodename).operator = 'normal';
end

% output to prmflow.image
prmflow.image = prmrecon2image(prmflow.image, prmflow.recon);

% imagehead fields claiming
BPfields = {'InstanceNumber', 'imagecenter', 'reconcenter', 'ShotNumber', 'SliceLocation', 'ImagePositionPatient'};
if ~isfield(prmflow.image, 'imageheadfields') || isempty(prmflow.image.imageheadfields)
    % set default image head fields
    prmflow.image.imageheadfields = BPfields;
else
    prmflow.image.imageheadfields = union(prmflow.image.imageheadfields, BPfields);
end

% % while with a preview
% if prmflow.recon.previewed
%     prmflow.image.imageheadfields = union(prmflow.image.imageheadfields, 'Sreconact');
% end

% check rebin
if isfield(prmflow, 'rebin')
    if ~isfield(prmflow.rebin, 'issloperebin') || ~prmflow.rebin.issloperebin
        if prmflow.protocol.gantrytilt~=0
            % It is a mistake!
            warning(['The reconstruction no longer support previous Axialrebin when gantry tilting! ' ...
                'Please replace the reconnode Axialrebin by the Sloperebin.']);
        end
    end
end

% private buffer
dataflow.buffer.(nodename) = struct();
% set initial ishot
if strcmpi(scan, 'axial')
    dataflow.buffer.(nodename).ishot = 0;
else
    dataflow.buffer.(nodename).ishot = 1;
end

% pipe line
if pipeline_onoff
    dataflow.pipepool.(nodename) = status.defaultpool;
    prmflow.pipe.(nodename).pipeline.kernellevel = 1;
    % max input
    if ~isavail(prmflow.pipe.(nodename).pipeline.inputmaxlimit) || (prmflow.pipe.(nodename).pipeline.inputmaxlimit < 1)
        prmflow.pipe.(nodename).pipeline.inputmaxlimit = 512;
    end
    % min output
    if ~isavail(prmflow.pipe.(nodename).pipeline.outputminlimit) || (prmflow.pipe.(nodename).pipeline.outputminlimit < 1)
        prmflow.pipe.(nodename).pipeline.outputminlimit = 1;
    end
    % max output
    if ~isavail(prmflow.pipe.(nodename).pipeline.outputmaxlimit) || (prmflow.pipe.(nodename).pipeline.outputmaxlimit < 1)
        inputlimit = prmflow.pipe.(nodename).pipeline.inputmaxlimit;
        switch scan
            % to make it tight
            case 'axial'
                prmflow.pipe.(nodename).pipeline.outputmaxlimit = ...
                    double(ceil(prmflow.recon.Nslice * prmflow.recon.delta_d / prmflow.recon.imagethickness) + 5); % int
            case {'helical', 'conveyor'}
                prmflow.pipe.(nodename).pipeline.outputmaxlimit = ...
                    double(ceil(inputlimit/prmflow.recon.Nviewprot*prmflow.recon.imagesperpitch) + 1 + ...
                    prmflow.recon.imagerely*2);     % int
            otherwise
                prmflow.pipe.(nodename).pipeline.outputmaxlimit = prmflow.system.defaultimagepoolsize;
        end
    end

    % clear the view counting
    prmflow.pipe.(nodename).pipeline.viewrescale = [0 1];
    % skip the normal priostep!
    prmflow.pipe.(nodename).pipeline.priostep = false;
    % We have a special priostep for BP

    % poolsize and poolsize ask
    if isavail(prmflow.pipe.(nodename).pipeline.currpoolsize) && prmflow.pipe.(nodename).pipeline.currpoolsize > 0
        currpoolsize = prmflow.pipe.(nodename).pipeline.currpoolsize;
    else
        currpoolsize = max(prmflow.system.defaultrawpoolsize, prmflow.pipe.(nodename).pipeline.inputmaxlimit);
    end
    if isavail(prmflow.pipe.(nodename).pipeline.nextpoolsize) && prmflow.pipe.(nodename).pipeline.nextpoolsize > 0
        nextpoolsize = prmflow.pipe.(nodename).pipeline.nextpoolsize;
    else
        nextpoolsize = max(prmflow.system.defaultimagepoolsize, prmflow.pipe.(nodename).pipeline.outputmaxlimit);
    end
    dataflow.pipepool.(nodename).poolsize = currpoolsize;
    prmflow.pipe.(nodename).pipeline.nextpoolsize = nextpoolsize;

    % ask objecttype and datasize for next node
    prmflow.pipe.(nodename).pipeline.nextobjecttype = 'image';
    prmflow.pipe.(nodename).pipeline.nextdatasize = double(prmflow.recon.imagesize(1) * prmflow.recon.imagesize(2));

    % default GPU on
    if prmflow.pipe.(nodename).pipeline.GPUonoff == -1
        prmflow.pipe.(nodename).pipeline.GPUonoff = 1;
    end

    if strcmpi(prmflow.recon.method, 'axial3d')
        % the nextshot is a pool to buffer half-slices rawdata
        dataflow.buffer.(nodename).nextshot = status.defaultpool;
        dataflow.buffer.(nodename).nextshot.poolsize = double(prmflow.recon.Nviewprot);
        dataflow.buffer.(nodename).nextshot.circulatemode = true;
        dataflow.buffer.(nodename).lastshotstart = true;
        dataflow.buffer.(nodename).lastshotend = false;
    else
        % for helical
        dataflow.buffer.(nodename).spaceshift = 0;
    end

    if strcmpi(prmflow.recon.scan, 'conveyor')
        % for conveyor, to set the viewexpand
        prmflow.pipe.(nodename).pipeline.viewexpand = prmflow.recon.viewexpand;
        % ini buffer
        dataflow.buffer.(nodename).Zview = [];
        dataflow.buffer.(nodename).Zsamp = [];
        dataflow.buffer.(nodename).Pz = 0;
    end

    if prmflow.recon.previewed
        % ask
        prmflow.pipe.(nodename).pipeline.nextobjecttype = 'sparseimage';
        % pool(2)
        nextnode = status.pipeline.(nodename).nextnode;
        dataflow.pipepool.(nodename)(2) = status.defaultpool;
        dataflow.pipepool.(nodename)(2).datafields = {'image', 'imagehead', 'Sreconact'};
        dataflow.pipepool.(nodename)(2).data = struct();
        dataflow.pipepool.(nodename)(2).poolsize = nextpoolsize;
        dataflow.pipepool.(nodename)(2).iscarried = true;
        dataflow.pipepool.(nodename)(2).carrynode = nextnode;
        dataflow.pipepool.(nodename)(2).carryindex = 1;
        % initial the pool data
        % (due to the dataflow.pipepool.(nodename)(2).iscarried is true, we need not to initial it.)
        if ~dataflow.pipepool.(nodename)(2).iscarried
            prevnode = status.pipeline.(nodename).prevnode;
            prevask = prmflow.pipe.(prevnode).pipeline;
            datasize = double(prmflow.recon.imagesize(1) * prmflow.recon.imagesize(2));
            databasis = prevask.nextdatabasis;
            poolsize = nextpoolsize;
            currgpuDevice = prevask.nextgpuDevice;
            dataflow.pipepool.(nodename)(2).data = initialpooldata(dataflow.pipepool.(nodename)(2).data, ...
                'sparseimage', poolsize, datasize, databasis, currgpuDevice);
        end
    end
else
    % ini image
    dataflow.image = zeros(prmflow.recon.imagesize(1)*prmflow.recon.imagesize(2), 0, 'single');
    dataflow.imagehead = struct();
end
% we will employ a setting in prmflow.recon to define the fields in imagehead, TBC.

% GPU on/off
prmflow = defaultGPUonoff(prmflow, status, nodename);
% while GPU on
if prmflow.pipe.(nodename).pipeline.GPUonoff > 0
    GPUfields = {'XY', 'upsampgamma', 'delta_d_up', 'midchannel_up', 'SID', 'ConeWeightScale', 'Zviewend', ...
        'delta_z', 'Zupsamp', 'filter', 'Zinterp', 'Zcrossinterp', 'forward', 'activeXY', ...
        'Spreview', 'previewkernel', 'previewcut'};
    prmflow.recon = putfieldsinGPU(prmflow.recon, GPUfields);
end

% % a test of blocked BP
% if strcmpi(scan, 'helical')
%     % prmflow.recon.XY = prmflow.recon.XY + 0.05;
%     prmflow.recon = woolrolls(prmflow.recon, prmflow.rebin, nodeprm);
% end

% CUDA on/off
if prmflow.pipe.(nodename).pipeline.CUDAonoff
    % call cuda
    if ismember(lower(prmflow.recon.method), {'helicalpiline', 'conveyorpiline', 'conveyor'})
        [dataflow, prmflow, status] = cudabackprojectionprepare(dataflow, prmflow, status);
    else
        warning('The CUDA BP function only support helical-piline yet.');
        prmflow.pipe.(nodename).pipeline.CUDAonoff = false;
    end
end

end

function cimage = prmrecon2image(cimage, recon)
% to output the image parameters after recon (to the post-recon nodes)
% call it after BP prepare.

cimage.Nimage = recon.Nimage;
cimage.imagesize = recon.imagesize;
cimage.center = recon.center;
cimage.imageincrement = recon.imageincrement;
cimage.imagethickness = recon.imagethickness;
cimage.reconmethod = recon.method;
cimage.voxelsize = recon.voxelsize;
cimage.delta_radius = recon.delta_d/recon.voxelsize;
cimage.effFOV = recon.effFOV;
cimage.databasis = recon.databasis;

% Note: the information of each image shall be read from imagehead.
end
