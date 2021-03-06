function [dataflow, prmflow, status] = reconnode_bonehardencali(dataflow, prmflow, status)
% recon node, bone-beamharden calibration
% [dataflow, prmflow, status] = reconnode_bonehardencali(dataflow, prmflow, status);

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
bonecaliprm = prmflow.pipe.(status.nodename);
% format version of calibration table
if isfield(bonecaliprm, 'corrversion')
    corrversion = bonecaliprm.corrversion;
else
    corrversion = 'v1.0';
end

% load beamharden
if isfield(bonecaliprm, 'beamhardencorr')
    % user defined a beamhardencorr file in pipe line
    beamhardencorr = loaddata(bonecaliprm.beamhardencorr, prmflow.IOstandard);
elseif isfield(dataflow, 'beamhardencorr')
    % it was done by a previous recon-node of beamhardencali (SUGGEST)
    beamhardencorr = dataflow.beamhardencorr;
elseif isfield(prmflow.corrtable, 'Beamharden')
    % it was load a Beamharden calibration table
    beamhardencorr = dataflow.Beamharden;
else
    error('Not found beamhardencorr to do bone-beamharden calibration!');
end

% to check is defined SYS
if ~isfield(prmflow, 'SYS')
    error('Not defined CT system in pipe line to bone-beamharden calibration!');
end
% NOTE: due to the reconnode_beamhardencali could change the prmflow.SYS.collimation.filter by adding an experiment based 
% effective filter, the reconnode_bonehardencali should be called after the reconnode_beamhardencali to employ that effective 
% filter.

% simu bone cali
bonehardencorr = simuBonehardencali(prmflow.SYS, beamhardencorr, corrversion);

% slice merge prm (I know the prmflow.SYS.detector has been merged)
bonehardencorr.startslice = prmflow.system.detector.startslice;
bonehardencorr.endslice = prmflow.system.detector.endslice;
bonehardencorr.mergescale = prmflow.system.detector.mergescale;
bonehardencorr.slicemerge = prmflow.system.detector.slicemerge;

% to return
dataflow.bonehardencorr = bonehardencorr{1};

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end