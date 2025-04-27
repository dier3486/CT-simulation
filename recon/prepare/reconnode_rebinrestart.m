function [dataflow, prmflow, status] = reconnode_rebinrestart(dataflow, prmflow, status)
% recon node, rebin restart
% [dataflow, prmflow, status] = reconnode_rebinrestart(dataflow, prmflow, status);

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

% left view number
prmflow.rebin.viewleft = prmflow.rebin.Nview - prmflow.rebin.viewread;

% update viewnumber and Nview
prmflow.rebin.viewnumber = prmflow.rebin.viewnumber - prmflow.rebin.viewread;
prmflow.rebin.Nview = min(prmflow.rebin.viewnumber, prmflow.rebin.maxviewnumber);

% output to recon
prmflow.recon.viewleft = prmflow.recon.Nview - prmflow.recon.viewread;
newNview = floor((prmflow.rebin.viewread + prmflow.rebin.Nview) / prmflow.rebin.Nfocal);
prmflow.recon.Nview = newNview - prmflow.recon.viewread;

% reset viewread
prmflow.rebin.viewread = 0;


% status
status.jobdone = true;

end