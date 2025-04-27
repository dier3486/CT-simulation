function [rawref, refblock] = airreference(rawdata, referrcut, refpixel, Nslice, flag_slice)
% air reference for rawdata, to use side module as refernce detector and 
% compare with referrcut
% [rawref, referr] = airreference2(rawdata, referrcut, refpixel, Npixel, Nslice, flag_slice);
% INPUT:
%   rawdata         raw data, in size Npixel*Nslice*Nview
%   referrcut       the cutoff by the air calibration
%   refpixel        number of reference pixels
%   Nslice          slice number
% OUPUT:
%   rawref          reference value (raw-air), in size 2*Nview
%   refblock        is blocked of the reference (referr>referrcut), in size 2*Nview
%                   which is used to judge if the refernce detectors are
%                   blocked by objects in FOV

% size or reshape
Nview = size(rawdata, 2);
rawdata = reshape(rawdata, [], Nslice, Nview);

% is slice independent
if nargin<5
    flag_slice = false;
end

% skip the edge slices
if Nslice>2 && ~flag_slice
    index_slice = 2:Nslice-1;
    Nref = Nslice-2;
else
    index_slice = 1:Nslice;
    Nref = Nslice;
end

% reference data
ref1 = reshape(rawdata(refpixel(1,:), index_slice, :), [], Nview);
ref2 = reshape(rawdata(refpixel(2,:), index_slice, :), [], Nview);

% reference
if ~flag_slice
    rawref = [mean(ref1, 1); mean(ref2, 1)];
else
    reflength = size(refpixel, 2);
    rawref = [squeeze(mean(reshape(ref1, reflength, Nref, Nview), 1)); ...
              squeeze(mean(reshape(ref2, reflength, Nref, Nview), 1))];
end

% reference error
referr = [std(ref1, 0, 1); std(ref2, 0 , 1)];

% % tmp debug
% referr = referr*22.8;

% reference blocked
refblock = false(size(referr));
Nfocal = size(referrcut, 2);
for ifocal = 1:Nfocal
    refblock(:, ifocal:Nfocal:end) = referr(:, ifocal:Nfocal:end) > max(referrcut(:, ifocal), 1e-5);
end

return

