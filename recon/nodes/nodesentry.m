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
        [dataflow, prmflow, status] = reconnode_beamhardencorr(dataflow, prmflow, status);
    case 'housefield'
        [dataflow, prmflow, status] = reconnode_housefieldcorr(dataflow, prmflow, status);
    case 'axialrebin'
        [dataflow, prmflow, status] = reconnode_Axialrebin(dataflow, prmflow, status);
    case 'filter'
        7;
    case 'backprojection'
        8;
    case 'statusmatrix'
        9;
    otherwise
        % handle
        myfun = str2func(['reconnode_' nodename]);
        [dataflow, prmflow, status] = myfun(dataflow, prmflow, status);
        % But we suggest to switch-case a node's name above as a register.
        % It will be easy to set a break for debug.
end

end