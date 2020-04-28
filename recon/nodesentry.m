function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename, varargin)
% call nodes
% [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename);

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
switch lower(nodename_slip{1})
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
    case 'housefield'
        % Housefield CT value correction
        [dataflow, prmflow, status] = reconnode_housefieldcorr(dataflow, prmflow, status);
    % rebin & FBP
    case 'axialrebin'
        % rebin for axial
        [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
    case 'filter'
        % filter
        [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);
    case {'backproject', 'backprojection', 'bp'}
        % back projection
        [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);
    case 'fbp'
        % temporary FBP function
        [dataflow, prmflow, status] = reconnode_FBPtmp(dataflow, prmflow, status);
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
    case 'nonelinearcali'
        % nonlinear calibration
        [dataflow, prmflow, status] = reconnode_nonelinearcali(dataflow, prmflow, status);
    case 'crosstalkcali'
        % crosstalk calibration
        [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status);
    case 'offfocalcali'
        % offfocal calibration
        [dataflow, prmflow, status] = reconnode_offfocalcali(dataflow, prmflow, status);
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

end