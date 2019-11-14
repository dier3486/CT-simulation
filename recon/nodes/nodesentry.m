function [dataflow, prmflow, status] = nodesentry(dataflow, prmflow, status, nodename)
% call nodes

switch nodename
    case 'initial'
        0;
    case 'loadrawdata'
        1;
    case 'loadcorrs'
        2;
    case 'log2'
        3;
    case 'aircorr'
        4;
    case 'HCcorr'
        5;
    case 'rebin'
        6;
    case 'filter'
        7;
    case 'backprojection'
        8;
    otherwise
        1;
end

end