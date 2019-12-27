function [rawref, referr] = airreference2(rawdata, airref, refpixel, Npixel, Nslice)
% air reference for rawdata, to use side module as refernce detector and 
% compare with the references in airref
% [rawref, referr] = airreference2(rawdata, airref, refpixel, Npixel, Nslice);
% INPUT:
%   rawdata         raw data, in size Npixel*Nslice*Nview
%   airref          air reference, in size (refpixel*2)*Nslice*Nview
%                   which has been interpolated from the sections to the
%                   views
%   refpixel        number of reference pixels
%   Npixel          pixel number of each slice
%   Nslice          slice number
% OUPUT:
%   rawref          reference value (raw-air), in size 2*Nview
%   referr          err of the reference, in size 2*Nview
%                   which is used to judge if the refernce detectors are
%                   blocked by objects in FOV

% size or reshape
if nargin<4
    [Npixel, Nslice, Nview] = size(rawdata);
else
    rawdata = reshape(rawdata, Npixel, Nslice, []);
    Nview = size(rawdata, 3);
end
airref = reshape(airref, 2*refpixel, Nslice, Nview);

% skip the edge slices
if Nslice>2
    index_slice = 2:Nslice-1;
    Nrefsl = Nslice-2;
else
    index_slice = 1:Nslice;
    Nrefsl = Nslice;
end

% referenece data
airref1 = airref(:, index_slice, :);
airref2 = airref(:, index_slice+Nslice, :);
ref1 = rawdata(1:refpixel, index_slice, :);
ref2 = rawdata(Npixel-refpixel+1:Npixel, index_slice, :);

% rawref
ref1 = -log2(mean(2.^reshape(-ref1, [], Nview), 1)) + log2(mean(2.^reshape(-airref1, [], Nview), 1));
ref2 = -log2(mean(2.^reshape(-ref2, [], Nview), 1)) + log2(mean(2.^reshape(-airref2, [], Nview), 1));
rawref = [ref1; ref2];

% reference error
ref1_err = reshape(ref1 - airref1, [], Nview);
ref1_err = reshape(ref1_err-mean(ref1_err, 1), Nrefsl, refpixel, Nview);
ref2_err = reshape(ref2 - airref2, [], Nview);
ref2_err = reshape(ref2_err-mean(ref2_err, 1), Nrefsl, refpixel, Nview);
% SVD 
referr = zeros(2, Nview);
for iview = 1:Nview
    s1 = svd(ref1_err(:, :, iview));
    referr(1, iview) = s1(1);
    s2 = svd(ref2_err(:, :, iview));
    referr(2, iview) = s2(1);
end

return

