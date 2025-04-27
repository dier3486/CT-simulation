function [airmain, referrcut, referenceKVmA] = ...
    aircalibration(rawdata, viewangle, refpixelindex, Nsection, Nfocal, KVmA)
% standard air calibration sub function v2
% [airmain, referrcut] = aircalibration(rawdata, viewangle, refpixelindex, Nsection, Nfocal);
% rawdata is raw data of air after log2, and reshaped by
% rawdata = reshape(rawdata, Npixel, Nslice, Nview);
% viewangle is the view angles, in [0, 2pi)
% refpixelindex is the index of the reference detectors on each side
% Nsection is the section number
% Nfocal is the focal spot number

% size
[Npixel, Nslice, Nview] = size(rawdata);
Np = Npixel*Nslice;

% reference
% skip the edge slices
if Nslice>2
    index_slice = 2:Nslice-1;
    Nrefslice = Nslice - 2;
else
    index_slice = 1:Nslice;
    Nrefslice = Nslice;
end
indexNp = reshape(1:Np, Npixel, Nslice);
refindex1 = indexNp(refpixelindex(1,:), index_slice);
refindex2 = indexNp(refpixelindex(2,:), index_slice);
airref1 = reshape(rawdata(refpixelindex(1,:), index_slice, :), [], Nview);
airref2 = reshape(rawdata(refpixelindex(2,:), index_slice, :), [], Nview);

% ini
airmain = zeros(Npixel*Nslice*Nfocal, Nsection);
referenceKVmA = zeros(Nfocal, Nsection);

% section angle range
delta = pi*2/Nsection;
sectangle = (0:Nsection).*delta;

% shift viewangle
viewangle = mod(viewangle + delta/2, pi*2);

% KVmA
if nargin > 5
    KVmA(KVmA<1.0) = 1.0;
    KVmA = -log2(KVmA);
else
    KVmA = zeros(1, Nview);
end

% loop the section
for isect = 1:Nsection
    % views in this section
    viewindex = (viewangle>=sectangle(isect)) & (viewangle<sectangle(isect+1));
    for ifocal = 1:Nfocal
        viewindex_ifocal = false(1, Nview);
        viewindex_ifocal(ifocal:Nfocal:end) = viewindex(ifocal:Nfocal:end);
        % airmain
        mainindex = (1:Np) + Np*(ifocal-1);
        airmain(mainindex, isect) = reshape(mean(rawdata(:, :, viewindex_ifocal), 3), [], 1);
        % air ref
        airref1(:, viewindex_ifocal) = airref1(:, viewindex_ifocal) - airmain(mainindex(refindex1), isect);
        airref2(:, viewindex_ifocal) = airref2(:, viewindex_ifocal) - airmain(mainindex(refindex2), isect);
        % mA ref
        referenceKVmA(ifocal, isect) = mean(KVmA(viewindex_ifocal));
    end
end

% std
refstd1 = std(airref1, 0, 1);
refstd2 = std(airref2, 0, 1);

% refernce error cut of block
referrcut = zeros(2, Nfocal);
for ifocal = 1:Nfocal
    referrcut(1, ifocal) = mean(refstd1(ifocal:Nfocal:end));
    referrcut(2, ifocal) = mean(refstd2(ifocal:Nfocal:end));
end

return