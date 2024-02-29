function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename, keepnode)
% call nodes
%   [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename);
% or, [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status);

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

if nargin<5
    keepnode = false;
end

% multi nodes
if nargin>=4 && iscell(nodename)
    for ii = length(nodename)
        [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename{ii}, keepnode);
    end
end

% set status.nodename
if nargin<4
    nodename = status.nodename;
elseif ~keepnode
    status.nodename = nodename;
    % if keepnode don't replace the status.nodename
end
% slip
nodename_slip = regexp(nodename, '_', 'split');

% switch nodes
try
    switch lower(nodename_slip{1})
        % IO and control
        case 'statusmatrix'
            % TBC
            0;
        case 'initial'
            % initial, what ever call this first
            [dataflow, prmflow, status] = reconinitial(dataflow, prmflow, status);
        case {'loadrawdata', 'readraw'}
            % read rawdata
            [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status);
        case 'loadcorrs'
            % load calibration tables
            [prmflow, status] = loadcalitables(prmflow, status);
        case 'pipelineprepare'
            % pipeline prepare
            [dataflow, prmflow, status] = reconnode_pipelineprepare(dataflow, prmflow, status);
            % pipeline prepare is not a preparenode :)
        % corrections
        case 'log2'
            % log2
            [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);
        case {'aircorr', 'air'}
            % air correction
            [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);
        case 'reference'
            % air reference correction
            [dataflow, prmflow, status] = reconnode_referencecorr(dataflow, prmflow, status);
        case 'badchannel'
            % fix badchannel
            [dataflow, prmflow, status] = reconnode_badchannelcorr(dataflow, prmflow, status);
        case {'crosstalk'}
            % crosstalk correction
            [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status);
        case {'beamharden', 'nonlinear'}
            % beam harden and nonlinear correction
            [dataflow, prmflow, status] = reconnode_nonlinearcorr(dataflow, prmflow, status);
        case 'offfocal'
            % off-focal correction
            [dataflow, prmflow, status] = reconnode_offfocalcorr(dataflow, prmflow, status);
        case 'fanangle'
            % interp to equal fanangle
            [dataflow, prmflow, status] = reconnode_fananglecorr(dataflow, prmflow, status);
        case {'hounsefield', 'housefield', 'hu'}
            % Hounsefield Units (HU) correction, (wrong spelling tolerenced for 'Housefield')
            [dataflow, prmflow, status] = reconnode_hounsefieldcorr(dataflow, prmflow, status);
        case {'materialdecomp', 'materialdecomposition', 'md'}
            % two-material decomposition
            [dataflow, prmflow, status] = reconnode_materialdecompcorr(dataflow, prmflow, status);
        % rebin & FBP
        case 'rowcombine'
            % row (slices) combine
            [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status);
        case 'axialrebin'
            % rebin for axial
            [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
            % plz use sloperebin in 3D recon
        case 'sloperebin'
            % new rebin method of axial
            [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);
        case 'helicalrebin'
            % Helical rebin
            [dataflow, prmflow, status] = reconnode_Helicalrebin(dataflow, prmflow, status);
        case 'upsample'
            % doule up sampling
            [dataflow, prmflow, status] = reconnode_Upsample(dataflow, prmflow, status);
        case 'filter'
            % filter
            [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);
        case {'backproject', 'backprojection', 'bp'}
            % back projection
            [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);
        case {'antiring', 'postantiring'}
            % anti ring artifact in image space
            [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);
        case 'boneharden'
            % bone harden correction 
            [dataflow, prmflow, status] = reconnode_bonehardencorr(dataflow, prmflow, status);
        case {'iteration', 'iterationrecon'}
            % iteration reconstrcution
            [dataflow, prmflow, status] = reconnode_Iterationrecon(dataflow, prmflow, status);
        % calibrations
        case 'aircali'
            % air calibration
            [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status);
        case 'beamhardencali'
            % beamharden calibration
            [dataflow, prmflow, status] = reconnode_beamhardencali(dataflow, prmflow, status);
        case 'bonehardencali'
            % bone-beamharden calibration
            [dataflow, prmflow, status] = reconnode_bonehardencali(dataflow, prmflow, status);
        case 'nonlinearcali'
            % nonlinear calibration
            [dataflow, prmflow, status] = reconnode_nonlinearcali(dataflow, prmflow, status);
        case 'crosstalkcali'
            % crosstalk calibration
            [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status);
        case 'offfocalcali'
            % offfocal calibration
            [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status);
        case 'detectorcali'
            % detector calibration
            [dataflow, prmflow, status] = reconnode_detectorcali(dataflow, prmflow, status);
        case 'ballcali'
            % Z-calibration
            [dataflow, prmflow, status] = reconnode_ballcali(dataflow, prmflow, status);
        case 'materialdecompcali'
            % (two)-material decompoistion calibration
            [dataflow, prmflow, status] = reconnode_materialdecompcali(dataflow, prmflow, status);
            % cali supports
        case 'inverserebin'
            % inverse the rebin, from parallel beams back to fan
            [dataflow, prmflow, status] = reconnode_inverserebin(dataflow, prmflow, status);
        case 'watergoback'
            % a calibration algorithm in getting ideal water projection
            [dataflow, prmflow, status] = reconnode_watergoback(dataflow, prmflow, status);
        case 'idealwater'
            % a calibration algorithm in getting ideal water projection
            [dataflow, prmflow, status] = reconnode_idealwater(dataflow, prmflow, status);
        case {'databackup', 'backup'}
            % backup inner data
            [dataflow, prmflow, status] = reconnode_databackup(dataflow, prmflow, status);
        case {'dataoutput', 'output'}
            % output images or calibration tables to file
            [dataflow, prmflow, status] = reconnode_dataoutput(dataflow, prmflow, status);
        case 'corrredirect'
            % copy dataflow.xxxcorr to prmflow.corrtable.xxx (used in online correction)
            [dataflow, prmflow, status] = reconnode_corrredirect(dataflow, prmflow, status);
        case 'datamean'
            % calculate the mean of rawdata
            [dataflow, prmflow, status] = reconnode_datamean(dataflow, prmflow, status);
        case 'systemconfigue'
            % call CT-simulation to configure system
            [dataflow, prmflow, status] = reconnode_Systemconfigure(dataflow, prmflow, status);
        case 'simulation'
            % call CT-simulation to do simulation
            [dataflow, prmflow, status] = reconnode_simulation(dataflow, prmflow, status);
        otherwise
            % function handle, call a function in name of reconnode_nodename
            myfunname = ['reconnode_' nodename_slip{1}];
            if any(exist(myfunname) == [2 5 6]) % can run
                myfun = str2func(myfunname);
                [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
            else
                status.jobdone = false;
                status.errorcode = -9;
                status.errormsg = sprintf('Not exist pipe line node function %s!', nodename);
            end
%             if regexp(nodename_slip{1}, 'prepare$')
%                 % a prepare node
%                 [prmflow, status] = myfun(prmflow, status);
%             else
%                 [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
%             end
            % It is a flexible way to include any recon nodes.
            % We suggest to register all the recon functions in the above cases to govern them.
            % Anyhow there are many 'hidden' nodes behind the explicit recon functions in running the pipeline, you may find
            % them in the folder ~/recon/supportnodes/.
    end
catch me
    status.jobdone = false;
    status.errorcode = -1;
    status.errormsg = me;
    if status.warning_onoff
        warning('error in calling pipe line node %s!', nodename);
    end
    % test
%     rethrow(me);
%     warning(me.getReport);
end

end