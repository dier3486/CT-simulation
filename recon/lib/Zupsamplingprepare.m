function interpstr = Zupsamplingprepare(recon, nodeprm, crossflag)
% to prepare the upsampling in Z interpolation

if nargin<3
    crossflag = 0;
end

% defualt coeff
if isfield(nodeprm, 'Gamma')
    Gamma = nodeprm.Gamma;
else
    Gamma = [1.0  0.0];
    % NOTE: Setting Gamma = [1.0  0.0] and Zupsampling = 1 will decay to linear interpolatopn.
end

if isfield(nodeprm, 'Zupsampling')
    interpstr.Zupsampling = nodeprm.Zupsampling;
else
    % default value is 1
    interpstr.Zupsampling = 1;
    % NOTE: A suggest Zupsampling value in Z-upsapmling is 4.
    % And, specially, to set Zupsampling=0 will omit the z-upsampling and repalce it by a continuous omiga4-interpolation.
end

if interpstr.Zupsampling >= 1
    % prepare a matrix for upsampling
    switch crossflag
        case 0
            % normal (for helical)
            interpstr.ZupMatrix = ntimesupmatrix(recon.Nslice, interpstr.Zupsampling, Gamma);
        case 1
            % cross shot, full
            A = ntimesupmatrix(recon.Nslice+4, interpstr.Zupsampling, Gamma);
            % 1st-shot
            interpstr.ZupMatrix_1 = A(3:end-2, :);
            interpstr.ZupMatrix_1(1, :) = interpstr.ZupMatrix_1(1, :) + A(1, :) + A(2, :);
            % end-shot
            interpstr.ZupMatrix_end = A(3:end-2, :);
            interpstr.ZupMatrix_end(end, :) = interpstr.ZupMatrix_1(end, :) + A(end, :) + A(end-1, :);
            % middle shots
            interpstr.ZupMatrix_mid = A(3:end-2, :);
        case 2
            % cross shot, half (for axial-iteration)
            % 0-shot
            A1 = ntimesupmatrix(recon.Nslice/2+2, interpstr.Zupsampling, Gamma);
            interpstr.ZupMatrix_0 = A1(2:end, :);
            interpstr.ZupMatrix_0(1, :) = interpstr.ZupMatrix_0(1, :) + A1(1, :);
            % end-shot
            interpstr.ZupMatrix_end = A1(1:end-1, :);
            interpstr.ZupMatrix_end(end, :) = interpstr.ZupMatrix_end(end, :) + A1(end, :);
            % middle shots
            interpstr.ZupMatrix_mid = ntimesupmatrix(recon.Nslice+2, interpstr.Zupsampling, Gamma);
        case 3
            % half pixel extrap
            A = ntimesupmatrix(recon.Nslice + 2, interpstr.Zupsampling, Gamma);
            A(2, :) = A(2, :) + A(1, :);
            A(end-1, :) = A(end-1, :) + A(end, :);
            m = floor(interpstr.Zupsampling/2);
            interpstr.ZupMatrix = A(2:end-1, m+1:end-m);
        otherwise
            % error
            1;
    end
    interpstr.Nfourp = 0;
else
    % prepare an 'omiga4table' for omiga4-interpolation to replace the z-upsampling
    if isfield(nodeprm, 'Zupsamptablesize') && ~isempty(nodeprm.Zupsamptablesize)
        Nfourp = nodeprm.Zupsamptablesize;
        interpstr = structmerge(interpstr, omiga4table(Gamma, Nfourp));
    else
        % the omiga4table is skipped, suppose to call the omiga4interp.m in z interpolation.
        interpstr.Nfourp = 0;
    end
end
interpstr.Cgamma = Gamma;

end