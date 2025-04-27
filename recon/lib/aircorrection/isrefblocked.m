function refblock = isrefblocked(refblock, refblk_new, blockwindow, scan, index_out, Do2i)
% if the reference blocked
% refblock = isrefblocked(refblock, refblk_new, index_out, blockwindow, Do2i, scan);
% It works like refblock = refblock | convolution(refblk_new, blockwindow).
% Though the alignment could be tough.

if nargin < 5 || isempty(index_out)
    index_out = 1 : size(refblock, 2);
end
if nargin < 6 || isempty(Do2i)
    Do2i = 0;
end

switch scan
    case 'axial'
        % convolution of refblk_new
        refblock_window = conv2(refblk_new, true(1, 2*blockwindow+1)) > 0;
        % refblock = refblock | refblock_window
        n = length(index_out);
        if n > blockwindow*2
            refblock(:, index_out(1: blockwindow*2)) = refblock(:, index_out(1: blockwindow*2)) | refblock_window(:, 1: blockwindow*2);
            refblock(:, index_out(blockwindow*2+1:end)) = refblock(:, index_out(blockwindow*2+1:end)) | refblock_window(:, blockwindow*2+1:end);
            % It's not an unnecessary move due to the index_out could be not unique.
        else
            refblock(:, index_out) = refblock(:, index_out) | refblock_window;
            % while the n is small enough the index_out is believed in unique.
        end
    case {'helical', 'halfaxial'}
        % convolution of refblk_new
        refblock_window = conv2(refblk_new, true(1, 2*blockwindow+1)) > 0;
        % to align the input-output index
        n = length(index_out);
        m = size(refblk_new, 2);
        [indexIn2, indexOut2] = indexinoutmap(m, n, [blockwindow, blockwindow], Do2i);
        % refblock = refblock | refblock_window
        refblock(:, index_out(indexOut2)) = refblock(:, index_out(indexOut2)) | refblock_window(:, indexIn2);
    case 'static'
        refblock(:, index_out) = refblock(:, index_out) | refblk_new;
    otherwise
        % was error
end

end