function detector = mergedetector(detector)
% slice merge a detector by itself
% detector = mergedetector(detector);
% warn: ireversible

% orig parametes in detector
Npixel = detector.Npixel;
Nslice = detector.Nslice;
slicemerge = detector.slicemerge;

% position
[detector.position, Nmergedslice] = detectorslicemerge(detector.position, Npixel, Nslice, slicemerge, 'mean');

% response
if isfield(detector, 'response') && size(detector.response, 1)>1
    [detector.response, ~] = detectorslicemerge(detector.response, Npixel, Nslice, slicemerge, 'mean');
end

% pixelarea
if isfield(detector, 'pixelarea') && size(detector.pixelarea, 1)>1
    [detector.pixelarea, ~] = detectorslicemerge(detector.pixelarea, Npixel, Nslice, slicemerge, 'sum');
end

% pixelnormal
if isfield(detector, 'pixelnormal') && size(detector.pixelnormal, 1)>1
    [detector.pixelnormal, ~] = detectorslicemerge(detector.pixelnormal, Npixel, Nslice, slicemerge, 'mean');
    % WARN: the mean of the normal vectors could be meaningless
end

% crossmatrix
if isfield(detector, 'crossmatrix')
    % TBC
    % are you sure to merge a sparse matrix?
    1;
end

% replace the parametes after merge
detector.Nslice = Nmergedslice;
detector.slicemerge = 1:Nmergedslice;

% But don't do this:
% detector.startslice = 1;
% detector.endslice = Nmergedslice;

return
