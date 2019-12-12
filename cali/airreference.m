function airref = airreference(rawair, refpixel, Npixel, Nslice)
% air reference, to use side module as refernce detector
% airref = airreference(rawair, refpixel, Npixel, Nslice);

% size or reshape
if nargin<3
    [Npixel, Nslice, Nsection] = size(rawair);
else
    rawair = reshape(rawair, Npixel, Nslice, []);
    Nsection = size(rawair, 3);
end
% skip the edge slices
if Nslice>2
    index_slice = 2:Nslice-1;
else
    index_slice = 1:Nslice;
end

ref1 = -log2(mean(2.^reshape(-rawair(1:refpixel, index_slice, :), [], Nsection), 1));
ref2 = -log2(mean(2.^reshape(-rawair(Npixel-refpixel+1:Npixel, index_slice, :), [], Nsection), 1));
airref = [ref1; ref2];

return

