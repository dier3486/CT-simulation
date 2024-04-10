function [dataflow, prmflow, status] = reconnode_Antiwindmill(dataflow, prmflow, status)
% recon node, anti windmill artifact on image
% TV based method
% [dataflow, prmflow, status] = reconnode_Antiwindmill(dataflow, prmflow, status);

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

% no images?
if ~isfield(dataflow, 'image')
    return;
end

% parameters from pipe
AWprm = prmflow.pipe.(status.nodename);

% GPU?
if isfield(status, 'GPUinfo') && ~isempty(status.GPUinfo)
    GPUonoff = true;
else
    GPUonoff = false;
end

% is real?
flag_real = isreal(dataflow.image);

% defualt parameter
[downsample, Gblur, TVmu, TVlambda, TVCrange, TVNiter, TVtol, TVlogC, fixlimit, fixsigma, boundaryOpt] = ...
    defaultprm(AWprm, flag_real);

% GPU Array
if GPUonoff
    [Gblur, TVmu, TVlambda, TVCrange, TVtol, TVlogC, fixlimit, fixsigma] = ...
        putinGPU(Gblur, TVmu, TVlambda, TVCrange, TVtol, TVlogC, fixlimit, fixsigma);
end

% image size
if isfield(prmflow, 'recon') && isfield(prmflow.recon, 'imagesize')
    imagesize = prmflow.recon.imagesize([2 1]);
else
    imagesize = size(dataflow.image, 1, 2);
end

% reshape
dataflow.image = reshape(dataflow.image, imagesize(1), imagesize(2), []);

if GPUonoff
    image0 = gpuArray(dataflow.image);
else
    image0 = dataflow.image;
end

% Gauss blur
if abs(Gblur) > eps
    image0 = gaussblur(dataflow.image, Gblur);
    % I know the gaussblur can treat complex image. 
end

% down sampling
if downsample > 1
    image0 = imagedown(image0, downsample);
end

% boundary option
switch lower(boundaryOpt)
    case {0, 'none'}
        0;
    case {1, 'copy'}
        image0 = cat(3, image0(:,:,1), image0, image0(:,:,end));
    otherwise
        1;
end

% expanded TV
DIM = 3;
if flag_real
    imageTV = GlogTV1D(image0, TVmu, TVlambda, TVlogC, [], DIM, TVCrange, TVNiter, TVtol);
else
    1;
    imageTV = GlogTV1D(real(image0), real(TVmu), real(TVlambda), real(TVlogC), [], DIM, real(TVCrange), ...
                real(TVNiter), real(TVtol)) + ...
              GlogTV1D(imag(image0), imag(TVmu), imag(TVlambda), imag(TVlogC), [], DIM, imag(TVCrange), ...
                imag(TVNiter), imag(TVtol)) .* 1i;
end

% diff
image0 = imageTV - image0;

% boundary option
switch lower(boundaryOpt)
    case {0, 'none'}
        0;
    case {1, 'copy'}
        image0 = image0(:,:,2:end-1);
    otherwise
        1;
end

% limit
% if length(fixlimit) > 1
%     sreal = real(image0) > 0;
%     simag = imag(image0) > 0;
%     fixlimit = abs(real(fixlimit(1))).*~sreal + abs(real(fixlimit(2))).*sreal + ...
%               (abs(imag(fixlimit(1))).*~simag + abs(imag(fixlimit(2))).*simag).*1i;
% end
if flag_real
    image0 = hardentanh(image0, fixsigma, fixlimit);
else
    image0 = hardentanh(real(image0), real(fixsigma), real(fixlimit)) + ...
             hardentanh(imag(image0), imag(fixsigma), imag(fixlimit)).*1i;
end

% up sampling
if downsample > 1
    image0 = imageupnGPU(image0, downsample);
    image0 = image0(1:imagesize(1), 1:imagesize(2), :);
end

% fix
dataflow.image = dataflow.image + gather(image0);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


% default parameters
function [downsample, Gblur, TVmu, TVlambda, TVCrange, TVNiter, TVtol, TVlogC, fixlimit, fixsigma, boundaryOpt] = ...
    defaultprm(AWprm, flag_real)

% down-sampling of the images in anti windmill
if isfield(AWprm, 'downsample')
    downsample = AWprm.downsample;
else
    % 2-times down-sampling
    downsample = 2;
end
% blur
if isfield(AWprm, 'Gblur')
    Gblur = AWprm.Gblur;
else
    Gblur = 0;
end
% TV 
if isfield(AWprm, 'TVmu')
    TVmu = AWprm.TVmu;
else
    TVmu = 0.03 + 0.03i.*~flag_real;
end
if isfield(AWprm, 'TVlambda')
    TVlambda = AWprm.TVlambda;
else
    TVlambda = TVmu./3;
end
if isfield(AWprm, 'TVCrange')
    TVCrange = AWprm.TVCrange;
else
    TVCrange = [-inf inf] + [-inf inf].*1i.*~flag_real;
end
if isfield(AWprm, 'TVNiter')
    TVNiter = AWprm.TVNiter;
else
    TVNiter = 100 + 100i.*~flag_real;
end
if isfield(AWprm, 'TVtol')
    TVtol = AWprm.TVtol;
else
    TVtol = 0.2 + 0.2i.*~flag_real;
end

if isfield(AWprm, 'TVlogC')
    TVlogC = AWprm.TVlogC;
else
    TVlogC = 4.0 + 4.0i.*~flag_real;
end
if isfield(AWprm, 'fixlimit')
    fixlimit = AWprm.fixlimit;
    % no inf plz
    if any(isinf(fixlimit))
        fixlimit(isinf(fixlimit)) = 1e9;
    end
else
    fixlimit = 100 + 100i.*~flag_real;
end
if isfield(AWprm, 'fixsigma')
    fixsigma = AWprm.fixsigma;
else
    fixsigma = 1e-3 + 1e-3i.*~flag_real;
end

if isfield(AWprm, 'boundaryOpt')
    boundaryOpt = AWprm.boundaryOpt;
else
    boundaryOpt = 'none';
end


end

