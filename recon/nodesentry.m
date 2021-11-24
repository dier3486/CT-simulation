function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename, varargin)
% call nodes
% [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename);

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

% to recurse
if nargin>4
    [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, [nodename varargin]);
end
if iscell(nodename)
    for ii = length(nodename)
        [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename{ii});
    end
end

% set status.nodename
status.nodename = nodename;
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
            [prmflow, status] = reconinitial(prmflow, status);
        case {'loadrawdata', 'readraw'}
            % read rawdata
            [dataflow, prmflow, status] = reconnode_readrawdata(dataflow, prmflow, status);
        case 'loadcorrs'
            % load calibration tables
            [prmflow, status] = loadcalitables(prmflow, status);
        % corrections
        case 'log2'
            % log2
            [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);
        case {'aircorr', 'air'}
            % air correction
            [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);
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
        % rebin & FBP
        case 'rowcombine'
            % row (slices) combine
            [dataflow, prmflow, status] = reconnode_Rowcombine(dataflow, prmflow, status);
        case 'axialrebin'
            % rebin for axial
            [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
        case 'sloperebin'
            % new rebin method of axial
            [dataflow, prmflow, status] = reconnode_Sloperebin(dataflow, prmflow, status);
        case 'filter'
            % filter
            [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);
        case {'backproject', 'backprojection', 'bp'}
            % back projection
            [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);
        case 'fbp'
            % temporary FBP function
            [dataflow, prmflow, status] = reconnode_FBPtmp(dataflow, prmflow, status);
        case {'antiring', 'postantiring'}
            % anti ring artifact in image space
            [dataflow, prmflow, status] = reconnode_Antiring(dataflow, prmflow, status);
        case 'boneharden'
            % bone harden correction 
            [dataflow, prmflow, status] = reconnode_bonehardencorr(dataflow, prmflow, status);
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
            myfun = str2func(['reconnode_' nodename_slip{1}]);
            [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
            % It is a flexible way to include any recon nodes.
            % But we suggest to register a node's name in above cases, that
            % will be easy to set breaks for debug.
    end
catch me
    status.jobdone = false;
    status.errorcode = -1;
    status.errormsg = me;
    warning('error in calling pipe line node %s', nodename);
    % test
%     rethrow(me);
%     warning(me.getReport);
end

end