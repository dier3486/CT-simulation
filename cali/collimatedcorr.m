function corrtable = collimatedcorr(corrtable, corrname, detector)
% collimator exposure of calibration tables
% corrtable = collimatedcorr(corrtable, corrname, detector);
% I know the detector has been slice merged

% max effective slice index
N = max(detector.endslice, corrtable.endslice);

% slice number
if isfield(corrtable, 'slicenumber')
    corrslicenum = corrtable.slicenumber;
else
    corrslicenum = corrtable.endslice - corrtable.startslice + 1;
end
% NOTE: corrslicenum is not always equal to corrtable.Nslice
if isfield(corrtable, 'slicemerge') && any(corrtable.slicemerge)
	corrslicemerge = corrtable.slicemerge;
elseif isfield(corrtable, 'mergescale') && corrtable.mergescale>0
	tmp = repmat(1 : corrslicenum/corrtable.mergescale, corrtable.mergescale, 1);
    corrslicemerge = tmp(:)';
else
    corrslicemerge = 1:corrslicenum;
end

% to mapping the index
indexbase = zeros(N, 1);
indexbase(corrtable.startslice : corrtable.endslice) = corrslicemerge;
indexbase = indexbase(detector.startslice : detector.endslice);

if isfield(corrtable, 'Nslice')
    Nslice = corrtable.Nslice;
else
    Nslice = corrslicenum/corrtable.mergescale;
end

slicemap = zeros(Nslice, 1);
slicemap(indexbase) = detector.slicemerge;

switch corrname
    case 'air'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'beamharden'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
        [corrtable.airrate, ~] = detectorslicemerge(corrtable.airrate, detector.Npixel, Nslice, slicemap, 'mean');
    case 'boneharden'
        1;
    case 'crosstalk'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'nonlinear'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'offfocal'
        [corrtable.airrate, Nmergedslice] = detectorslicemerge(corrtable.airrate, detector.Npixel, Nslice, slicemap, 'mean');
    case 'idealwater'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    otherwise
        % do nothing
        1;
end
corrtable.Nslice = Nmergedslice;

corrtable.startslice = detector.startslice;
corrtable.endslice = detector.endslice;
corrtable.mergescale = detector.mergescale;
corrtable.slicemerge = detector.slicemerge;
    
end