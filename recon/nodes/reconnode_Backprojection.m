function [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status)
% recon node, BP (after filter)
% [dataflow, prmflow, status] = reconnode_Backprojection(dataflow, prmflow, status);

% BP prepare
[dataflow, prmflow, status] = reconnode_BPprepare(dataflow, prmflow, status);

% I know the prmflow.recon.method was filled by reconnode_BPprepare
switch lower(prmflow.recon.method)
    case 'axial2d'
        % 2D Axial
        [dataflow, prmflow, status] = reconnode_Axial2DBackprojection(dataflow, prmflow, status);
    case 'axial3d'
        % 3D Axial
        [dataflow, prmflow, status] = reconnode_Axial3DBackprojection(dataflow, prmflow, status);
    otherwise
        error('Sorry, the reconstruction %s is not supported yet!', prmflow.recon.method);
end

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end