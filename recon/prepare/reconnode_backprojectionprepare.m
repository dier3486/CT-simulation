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

% common prepare
prmflow.recon = commonbpprepare(prmflow.recon, prmflow.protocol, prmflow.system, nodeprm);

% scan
scan = lower(prmflow.recon.scan);
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
        otherwise
            recon_method = '';
    end
end

if ~strncmpi(recon_method, prmflow.recon.scan, length(prmflow.recon.scan))
    prmflow.recon.method = [prmflow.recon.scan recon_method];
else
    prmflow.recon.method = recon_method;
end
% The default recon method of axial is Axial2D

% switch recon method
switch lower(prmflow.recon.method)
    case {'axial2d', 'axial3d'}
        prmflow.recon = axial3Dprepare(prmflow.recon, nodeprm);
    case {'helical', 'helical3d', 'helicalpiline'}
        % helical is always 3D
        prmflow.recon = helicalprepare(prmflow.recon, nodeprm);
    case {'axialhalf'}
        0;
        % TBC
    case {'cradle', 'cradlepiline'}
        % cradle(piline)
        prmflow.recon = cradleprepare(prmflow.recon, nodeprm);
    otherwise
        % do nothing
        1;
        % no topo
end

% to single
prmflow.recon = everything2single(prmflow.recon, 'double', 'single');

% output to prmflow.image
prmflow.image = prmrecon2image(prmflow.image, prmflow.recon);

% GPU
GPUfields = {'XY', 'upsampgamma', 'delta_d_up', 'midchannel_up', 'SID', 'ConeWeightScale', 'Zviewend', ...
    'imageincrement', 'delta_z', 'imageZgrid', 'Zupsamp', 'filter', 'Zinterp', 'Zcrossinterp', 'forward'};
prmflow.recon = putfieldsinGPU(prmflow.recon, GPUfields);

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
    % minoutput
    if ~isfield(prmflow.pipe.(nodename).pipeline, 'outputminlimit')
        prmflow.pipe.(nodename).pipeline.outputminlimit = 1;
    end
    % skip the normal priostep!
    prmflow.pipe.(nodename).pipeline.priostep = false;
    % We have a special priostep for BP

    % ask objecttype and datasize for next node
    prmflow.pipe.(nodename).pipeline.nextobjecttype = 'image';
    prmflow.pipe.(nodename).pipeline.nextdatasize = double(prmflow.recon.imagesize(1) * prmflow.recon.imagesize(2));
    
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
else
    % ini image
    dataflow.image = zeros(prmflow.recon.imagesize(1)*prmflow.recon.imagesize(2), 0, 'single');
    dataflow.imagehead = struct();
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
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
% Note: the information of each image shall be read from imagehead.

end
