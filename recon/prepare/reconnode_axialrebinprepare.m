function [dataflow, prmflow, status] = reconnode_axialrebinprepare(dataflow, prmflow, status)
% prepare node, axial-2D rebin prepare
% [dataflow, prmflow, status] = reconnode_axialrebinprepare(dataflow, prmflow, status);

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

% sure?
prmflow.raw.scan = 'axial';

% no slope
nodename = status.nodename;
prmflow.pipe.(nodename).sloperebin = false;

% rebin prepare
[dataflow, prmflow, status] = reconnode_rebinprepare(dataflow, prmflow, status);

end