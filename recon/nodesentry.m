function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename)
% call nodes

status.nodename = nodename;
switch lower(nodename)
    case 'initial'
        % initial, what ever call this first
        [prmflow, status] = reconinitial(status);
    case {'loadrawdata', 'readraw'}
        % read rawdata
        [dataflow, prmflow, status] = readrawdata(status.reconcfg, dataflow, prmflow, status);
    case 'loadcorrs'
        % load calibration tables
        [prmflow, status] = loadcalitables(prmflow, status);
    case 'log2'
        % log2
        [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);
    case {'aircorr', 'air'}
        % air correction
        [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);
    case 'beamharden'
        % beam harden correction
        [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status);
    case 'housefield'
        % Housefield CT value correction
        [dataflow, prmflow, status] = reconnode_housefieldcorr(dataflow, prmflow, status);
    case 'axialrebin'
        % rebin for axial
        [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
    case 'filter'
        % filter
        [dataflow, prmflow, status] = reconnode_Filter(dataflow, prmflow, status);
    case 'backprojection'
        8;
    case 'fbp'
        % temporary FBP function
        [dataflow, prmflow, status] = reconnode_FBPtmp(dataflow, prmflow, status);
    case 'statusmatrix'
        9;
    case 'aircali'
        % air calibration
        [dataflow, prmflow, status] = reconnode_aircali(dataflow, prmflow, status);
    otherwise
        % function handle, call a function in name of reconnode_nodename
        myfun = str2func(['reconnode_' nodename]);
        [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
        % It is a flexible way to include any recon nodes.
        % But we suggest to register a node's name in above cases, that 
        % will be easy to set breaks for debug.
end

end