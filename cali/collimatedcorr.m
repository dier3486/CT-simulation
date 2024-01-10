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

if isfield(corrtable, 'Nslice')
    Nslice = corrtable.Nslice;
else
    Nslice = corrslicenum/corrtable.mergescale;
end

if Nslice == detector.Nslice
    % simple case
    slicemap = (1:Nslice)';
elseif Nslice == 1
    % slice independent
    slicemap = 1;
else
    % to mapping the index
    indexbase = zeros(N, 1);
    indexbase(corrtable.startslice : corrtable.endslice) = corrslicemerge;
    indexbase = indexbase(detector.startslice : detector.endslice);
    slicemap = zeros(Nslice, 1);
    slicemap(indexbase) = detector.slicemerge;
end

Nmergedslice = Nslice;
switch corrname
    case 'air'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'beamharden'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
        if isfield(corrtable, 'airrate')
            [corrtable.airrate, ~] = detectorslicemerge(corrtable.airrate, detector.Npixel, Nslice, slicemap, 'mean');
        end
    case {'boneharden', 'iteration'}
        if Nslice > 1
            [corrtable.effbeamfilter, Nmergedslice] = detectorslicemerge(corrtable.effbeamfilter, detector.Npixel, Nslice, slicemap, 'mean');
            [corrtable.curvematrix, ~] = detectorslicemerge(corrtable.curvematrix, 1, Nslice, slicemap, 'mean');
        else
            % I know in slice independent method of bone-harden, the Nslice is 1
            Nmergedslice = Nslice;
        end
    case 'crosstalk'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'nonlinear'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
    case 'offfocal'
        if isfield(corrtable, 'airrate')
            [corrtable.airrate, Nmergedslice] = detectorslicemerge(corrtable.airrate, detector.Npixel, Nslice, slicemap, 'mean');
        end
    case 'idealwater'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, detector.Npixel, Nslice, slicemap, 'mean');
        [corrtable.indexrange, ~] = detectorslicemerge(corrtable.main, 2, Nslice, slicemap, 'roundmean');
    case 'hounsefield'
        [corrtable.main, Nmergedslice] = detectorslicemerge(corrtable.main, 1, Nslice, slicemap, 'mean');
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