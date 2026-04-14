function [dataflow, prmflow, status] = reconnode_smartaircali(dataflow, prmflow, status)
% air calibration
% [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)

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
Npixel = prmflow.raw.Npixel;
Nslice = prmflow.raw.Nslice;
Nfocal = prmflow.raw.Nfocal;
Nview = prmflow.raw.Nview;
Nshot = prmflow.raw.Nshot;
Nviewpershot = prmflow.raw.viewpershot(1);
Nviewprot = prmflow.raw.Nviewprot;
Nmulti = Nviewpershot/Nviewprot;
angulationcode = prmflow.system.angulationcode;
angulationzero = prmflow.system.angulationzero;     % double
anglepercode = (pi*2) / double(angulationcode);

% parameters to use
caliprm = prmflow.pipe.(status.nodename);
if isfield(caliprm, 'Nsection')
    Nsection = caliprm.Nsection;
else
    Nsection = 1;
end
if isfield(caliprm, 'refpixel')
    refpixel = caliprm.refpixel;
else
    refpixel = 0;
end
if isfield(caliprm, 'refpixelskip')
    refpixelskip = caliprm.refpixelskip;
else
    refpixelskip = 0;
end
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    % default version is v1.0
    corrversion = 'v1.0';
end
if isfield(caliprm, 'referrcutscale')
    referrcutscale = caliprm.referrcutscale;
else
    referrcutscale = 1;
end
if isfield(caliprm, 'stabletol')
    stabletol = caliprm.stabletol;
else
    stabletol = 0.05;
end

% paramters for smartaircorr
smartaircorr = caliprmforcorr(prmflow, corrversion);
smartaircorr.Nslice = Nslice;
smartaircorr.Nsection = Nsection;
smartaircorr.Nview = Nview;
smartaircorr.refpixel = refpixel;
smartaircorr.refnumber = 0;

% viewangle
viewangle = double(dataflow.rawhead.AngleEncoder).*anglepercode + angulationzero;

% rawdata mean of multi shots
dataflow.rawdata = mean(reshape(dataflow.rawdata, Npixel, Nslice, Nviewpershot, Nshot), 4);
% KVmA
KVmA = dataflow.rawhead.mA.*dataflow.rawhead.KV;
% check if the raw is stable
if stablecheck(dataflow.rawdata, stabletol)
    status.errormsg = 'The air calibration is banned due to the unstable air data! please redo the data scan.';
%     warning(status.errormsg);
    status.jobdone = true;
    status.errorcode = 2;
    return;
end

dataflow.rawdata = reshape(dataflow.rawdata, [], Nview);
[~, index_min] = min(viewangle);
dataflow.rawdata = [dataflow.rawdata(:,index_min:end),dataflow.rawdata(:,1:index_min-1)];
if isfield(prmflow.corrtable,'smartaircali') || isfield(prmflow.corrtable,'Smartaircali')
    airraw = reshape(prmflow.corrtable.Smartaircali.airraw, [], Nview);
    airmask = reshape(prmflow.corrtable.Smartaircali.airmask, [], Nview);
    smartaircorr.main = sum((dataflow.rawdata - airraw).*airmask, 2)./sum(airmask,2);
    smartaircorr.airraw = prmflow.corrtable.Smartaircali.airraw;
    smartaircorr.airmask = prmflow.corrtable.Smartaircali.airmask;
else
    % to create the original smartair correction table using the smartair
    % rawdata acquired on the same day as the normal air data
    smartaircorr.main = zeros(Npixel, Nslice, 'single');
    smartaircorr.airraw = dataflow.rawdata;
    w = weight(dataflow.rawdata);
    smartaircorr.airmask = w;
end

% bad channel warning
index_badchannel = find(abs(smartaircorr.main)>0.2);
if size(index_badchannel) > 0
    fprintf("Bad channel index is");
    fprintf("%d ", index_badchannel);
end

% mainsize
smartaircorr.mainsize = length(smartaircorr.main(:));

% to return
dataflow.smartaircorr = smartaircorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function r = stablecheck(rawdata, tol)

rawmean = 2.^mean(reshape(rawdata, [], length(rawdata)));
r = any( abs(rawmean./mean(rawmean) - 1) > tol);

end


function w = weight(rawdata)
    w = ones(size(rawdata), 'logical');
    dif0 = diff(rawdata, 1, 2);
    mean0 = mean(abs(dif0),2);
    dif1 = [zeros(size(rawdata,1), 1), dif0];
    dif2 = [dif0, zeros(size(rawdata,1), 1)];
    w(abs(dif1) > mean0*5 ) = 0;
    w(abs(dif2) > mean0*5 ) = 0;
    w(rawdata>-1.5) = 0;
end
