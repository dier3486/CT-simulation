function w = weightofslicemerge(detector)
% weighting in slice merging
% w = weightofslicemerge(detector);

w = detectorslicemerge(ones(detector.Nslice, 1), 1, detector.Nslice, detector.slicemerge, 'sum');
w = w(detector.slicemerge);
w = reshape(1./w.*detector.mergescale, 1, detector.Nslice);

end