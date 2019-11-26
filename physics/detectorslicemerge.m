function [Pout, Nmergedslice] = detectorslicemerge(Pin, detector, mergemethod)
% merge the slices

if nargin<3
    mergemethod = 'sum';
end

Nslice = detector.Nslice;
Npixel = detector.Npixel;

if all(detector.slicemerge == 1:Nslice)
    % skip
    Pout = Pin;
    Nmergedslice = Nslice;
    return;
end

% to merge
Nmergedslice = max(detector.slicemerge);
Pin = reshape(Pin, Npixel, Nslice, []);
Pout = zeros(Npixel, Nmergedslice, size(Pin, 3));
mergeweight = zero(1, Nmergedslice);
for ii = 1:Nslice
    index_ii = detector.slicemerge(ii);
    Pout(:, index_ii, :) = Pout(:, index_ii, :) + Pin(:, ii, :);
    mergeweight(index_ii) = mergeweight(index_ii) + 1;
end
switch lower(mergemethod)
    case 'mean'
        Pout = Pout./mergeweight;
    case 'sum'
        % do nothing
        1;
    otherwise
        % sum
        1;
end
Pout = reshape(Pout, Npixel*Nmergedslice, []);

end