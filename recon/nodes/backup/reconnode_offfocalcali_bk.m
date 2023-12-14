function [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status)
% recon node, off-focal calibration
% [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status);
% NOTE: this code only did some copy/paste works, the real off-focal calibration is manually adjusting the parameters, 
% Good Luck.

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
caliprm = prmflow.pipe.(status.nodename);

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end

% off-focal kernel parameter
if isfield(caliprm, 'offfocalkernel')
    offfocalkernel = caliprm.offfocalkernel;
elseif isfield(prmflow.system, 'offfocalkernel')
    offfocalkernel = prmflow.system.offfocalkernel;
else
    offfocalkernel = [];
end
% I know the offfocalkernel .xml shall be defined in reconxml.system or reconxml.pipe.offfocalcali

% debug
% offfocalkernel = myxml2struct('E:\PANGU.DAT\CT\SINO\Config\bhcali\offfocalparameter.xml');

% load offfocalkernel file
if ischar(offfocalkernel)  && ~isempty(offfocalkernel)
    offfocalkernel = readcfgfile(offfocalkernel);
end

% find out the coupled offfocal kernel
offcorrprm = offfocalloadkernel(offfocalkernel, prmflow.protocol);
% the codes are:
% offcorrprm = struct();
% if ~isempty(offfocalkernel)
%     % cell to list
%     offfocalkernel = structcellcurse(offfocalkernel);
%     % check KV & bowtie
%     isKV = [offfocalkernel.offfocalpara(:).KVp] == prmflow.protocol.KV;
%     isbowtie = strcmpi({offfocalkernel.offfocalpara(:).Bowtietype}, prmflow.protocol.bowtie);
%     % found any?
%     index_cp = find(isKV & isbowtie, 1);
%     if ~isempty(index_cp)
%         % check collimator
%         iscollimator = strcmp({offfocalkernel.offfocalpara(index_cp).collimation(:).collimationwidth}, ...
%             prmflow.protocol.collimator);
%         index_colli = find(iscollimator, 1);
%         if ~isempty(index_colli)
%             % found
%             offcorrprm = offfocalkernel.offfocalpara(index_cp).collimation(index_colli);
%         end
%     end
% end

% to check if we load an offfocal corr table
if isfield(prmflow.corrtable, 'Offfocal')
    offfocalbase = prmflow.corrtable.Offfocal;
elseif isfield(prmflow.corrtable, 'Beamharden')
    offfocalbase = prmflow.corrtable.Beamharden;
elseif isfield(dataflow, 'beamhardencorr')
    offfocalbase = dataflow.beamhardencorr;
else
    error('No off-focal base line found in prmflow and dataflow to do calibration!');
end

% offfocalcorr
offfocalcorr = caliprmforcorr(prmflow, corrversion);
% merge
offfocalcorr = structmerge(offfocalcorr, offcorrprm);
offfocalcorr = structmerge(offfocalcorr, offfocalbase);

% steps
if isfield(offfocalcorr, 'offintensity')
    offfocalcorr.stepnumber = length(offfocalcorr.offintensity);
else
    offfocalcorr.stepnumber = 1;
end

% to return 
dataflow.offfocalcorr =  offfocalcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end