function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename)
% call nodes

status.nodename = nodename;
switch lower(nodename)
    case 'initial'
        [prmflow, status] = reconinitial(status);
    case {'loadrawdata', 'readraw'}
        [dataflow, prmflow, status] = readrawdata(status.reconcfg, dataflow, prmflow, status);
    case 'loadcorrs'
        [prmflow, status] = loadcalitables(prmflow, status);
    case 'log2'
        [dataflow, prmflow, status] = reconnode_log2(dataflow, prmflow, status);
    case {'aircorr', 'air'}
        [dataflow, prmflow, status] = reconnode_aircorr(dataflow, prmflow, status);
    case 'hccorr'
        5;
    case 'rebin'
        6;
    case 'filter'
        7;
    case 'backprojection'
        8;
    case 'statusmatrix'
        9;
    otherwise
        1;
end

end