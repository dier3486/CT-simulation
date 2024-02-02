function [dataflow, prmflow, status] = reconnode_materialzmap(dataflow, prmflow, status)
% recon node, after two-material decomposed reconstruction create the Z-map
% [dataflow, prmflow, status] = reconnode_materialzmap(dataflow, prmflow, status);

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

% parameters to use in prmflow
Nimage = prmflow.recon.Nimage;
imagesize = prmflow.recon.imagesize;

% I know the material decomposition calibration table is
mdcorr = prmflow.corrtable.Materialdecomp;
% tmply used, which shall be moved to prepare (of reconnode_materialdecompcorr).
ZA = mdcorr.ZA;
ZB = mdcorr.ZB;
% tmply hard coded paramters
rhoaircut = 20;
Rneighbor = 10;
sigmaneib = Rneighbor*2;
sigmarho = 0.1;

rhomap = reshape(real(dataflow.image), [], Nimage);
Rmap = reshape(imag(dataflow.image)./rhomap, [], Nimage);
Zmap0 = Ramp.*(ZB^3 - ZA^3) + ZA^3;
Zmap0(rhomap < rhoaircut) = 0;
Zmap0 = abs(Zmap0).^(1/3).*sign(Zmap0);

[X, Y] = meshgrid(1:imagesize(1), 1:imagesize(2));
X = X(:); Y = Y(:);
for ii = 1 : prod(imagesize)
    r2 = (X - X(ii)).^2 + (Y - Y(ii)).^2;
    s = sqrt(r2) <= Rneighbor;
    
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end
