function w = weightofslicemerge(detector)
% weighting in slice merging
% w = weightofslicemerge(detector);

w = detectorslicemerge(ones(1, detector.Nslice), 1, detector.Nslice, detector.slicemerge, 'sum');
w = w(detector.slicemerge);
w = 1./w.*detector.mergescale;

end