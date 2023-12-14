function refblock = isrefblocked(refblock, prmflow, referr, startviewindex)
% if the reference blocked

if nargin<4
    startviewindex = [];
end

blockwindow = prmflow.raw.air.blockwindow;
referrcut = prmflow.raw.air.referrcut;
Nviewprot = prmflow.raw.Nviewprot;
% Nrefblock = size(refblock, 2);

Nview_ref = size(referr, 2);
refblock0 = referr > max(referrcut(:), 1e-5);
refblock_window = conv2(refblock0, true(1, 2*blockwindow+1)) > 0;


if isempty(startviewindex) || isempty(refblock)
    % single data block
    switch prmflow.raw.scan
        case 'static'
            refblock = refblock0;
        case 'helical'
            refblock = refblock_window(:, blockwindow+1:end-blockwindow);
        case 'axial'
            refblock = refblock_window(:, blockwindow+1:end-blockwindow);
            % periodic boundary for per shot
            refblock(:, 1 : blockwindow) = ...
                refblock(:, 1 : blockwindow) | refblock_window(:, end-blockwindow+1 : end);
            refblock(:, end-blockwindow+1 : end) = ...
                refblock(:, end-blockwindow+1 : end) | refblock_window(:, 1 : blockwindow);
        otherwise
            error('Can not run the air correction for %s scan!', prmflow.raw.scan);
    end
else
    % multi data blocks (pipe-line replication)
    % Note: in this mode, the refblock must be initialized.
    switch prmflow.raw.scan
        case 'static'
            refblock(:, startviewindex:startviewindex+Nview_ref-1) = refblock0;
        case 'helical'
            % cum up
            refblock(:, startviewindex:startviewindex+Nview_ref+blockwindow-1) = ...
                refblock(:, startviewindex:startviewindex+Nview_ref+blockwindow-1) | refblock_window(:, blockwindow+1:end);
        case 'axial'
            % periodic boundary
            index = mod((startviewindex-blockwindow:startviewindex+Nview_ref+blockwindow-1) - 1, Nviewprot) + 1;
            refblock(:, index(blockwindow+1 : end-blockwindow)) = refblock(:, index(blockwindow+1 : end-blockwindow)) ...
                | refblock_window(:, blockwindow+1:end-blockwindow);
            refblock(:, index(1:blockwindow)) = refblock(:, index(1:blockwindow)) ...
                | refblock_window(:, 1:blockwindow);
            refblock(:, index(end-blockwindow+1:end)) = refblock(:, index(end-blockwindow+1:end)) ...
                | refblock_window(:, end-blockwindow+1:end);
            % Note, it seems equal to refblock(:, index) |= refblock_window, but exactly not.
        otherwise
            error('Can not run the air correction for %s scan!', prmflow.raw.scan);
    end

end


end