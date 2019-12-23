function [airmain, airref] = aircalibration(rawdata, viewangle, refpixel, Nsection, Nfocal)
% standard air calibration sub function
% [airmain, airref] = aircalibration(rawdata, viewangle, refpixel, Nsection, Nfocal);
% rawdata is raw data of air after log2, and reshaped by
% rawdata = reshape(rawdata, Npixel, Nslice, Nview);
% viewangle is the view angles, in [0, 2pi)
% refpixel is the number of the reference detectors on each side
% Nsection is the section number
% Nfocal is the focal spot number

% size
[Npixel, Nslice, Nview] = size(rawdata);
Np = Npixel*Nslice;

% ini
airmain = zeros(Npixel*Nslice*Nfocal, Nsection);
airref = zeros(2*Nfocal, Nsection);

% section angle range
delta = pi*2/Nsection;
sectangle = (0:Nsection).*delta;

% shift viewangle
viewangle = mod(viewangle + delta/2, pi*2); 

% airref of each view
airref_all = airreference(rawdata, refpixel);

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
        % airref
        refindex = (1:2) + 2*(ifocal-1);
        airref(refindex, isect) = mean(airref_all(:, viewindex_ifocal), 2);
    end
end

return