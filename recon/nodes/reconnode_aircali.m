function [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status)
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
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nfocal = prmflow.recon.Nfocal;
Nview = prmflow.recon.Nview;

% parameters to use
caliprm = prmflow.pipe.(status.nodename);
if isfield(caliprm, 'Nsection')
    Nsection = caliprm.Nsection;
else
    Nsection = 24;
end
if isfield(caliprm, 'refpixel')
    refpixel = caliprm.refpixel;
else
    refpixel = 16;
end
if isfield(caliprm, 'firstangle')
    firstangle = caliprm.firstangle;
else
    firstangle = 0;
end
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    % default version is v1.0
    corrversion = 'v1.0';
end
if isfield(caliprm, 'stabletol')
    stabletol = caliprm.stabletol;
else
    stabletol = 0.05;
end

% paramters for aircorr
aircorr = caliprmforcorr(prmflow, corrversion);
aircorr.Nslice = Nslice;
aircorr.Nsection = Nsection;
aircorr.firstangle = firstangle;
aircorr.refpixel = refpixel;
aircorr.refnumber = 2;

% shift viewangle
viewangle = dataflow.rawhead.viewangle - firstangle;
% reshape and KVmA
dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, Nview);
KVmA = dataflow.rawhead.mA.*dataflow.rawhead.KV;
% check if the raw is stable
if stablecheck(dataflow.rawdata, stabletol)
    status.errormsg = 'The air calibration is banned due to the unstable air data! please redo the data scan.';
%     warning(status.errormsg);
    status.jobdone = true;
    status.errorcode = 2;
    return;
end
% airmain and reference
% v1
% [aircorr.main, aircorr.reference] = aircalibration(dataflow.rawdata, viewangle, refpixel, Nsection, Nfocal);
% v2
[aircorr.main, aircorr.referrcut, aircorr.referenceKVmA] = ...
    aircalibration2(dataflow.rawdata, viewangle, refpixel, Nsection, Nfocal, KVmA);
% mainsize
aircorr.mainsize = length(aircorr.main(:));

% to return
dataflow.aircorr = aircorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end


function r = stablecheck(rawdata, tol)

rawmean = 2.^mean(reshape(rawdata, [], length(rawdata)));
r = any( abs(rawmean./mean(rawmean) - 1) > tol);

end
