function detector = collimatorexposure(collimator, detector, det_corr, collimatorexplain)
% explain the collimator

% prepare
det_corr.position = reshape(det_corr.position, det_corr.Npixel, det_corr.Nslice, 3);
detector.Npixel = det_corr.Npixel;

% no collimatorexplain?
if nargin<4
    collimatorexplain = [];
end
if isempty(collimatorexplain)
    % load default explain (hard code)
    detector = hardcodeexposure(collimator, detector, det_corr);
    return
end

% explain the collimator
Nc = length(collimatorexplain(:));
for icoll = 1:Nc
    coll = collimatorexplain{icoll};
    if strcmpi(collimator, coll.name)
        % get the collimator protocol
        Nslice_orig = size(detector.position, 2);
        detector.Nslice = coll.Nslice;
        if isfield(coll, 'startslice') && isfield(coll, 'endslice')
            detector.startslice = coll.startslice;
            detector.endslice = coll.endslice;
        else
            detector.startslice = (Nslice_orig - coll.Nslice)/2 + 1;
            detector.endslice = (Nslice_orig + coll.Nslice)/2;
        end
        if isfield(coll, 'slicemerge')
            detector.slicemerge = coll.slicemerge;
        elseif isfield(coll, 'mergescale')
            tmp = repmat(1:coll.Nslice/coll.mergescale, coll.mergescale, 1);
            detector.slicemerge = tmp(:)';
        else
            detector.slicemerge = 1:detector.Nslice;
        end
        if isfield(coll,  'mergescale')
            detector.mergescale = coll.mergescale;
        else
            detector.mergescale = 1;
        end
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO * detector.mergescale;
        return
    end
end

% why go to here?
if strcmpi(collimator, 'all') || strcmpi(collimator, 'open')
    % all on for lazy
    detector.position = reshape(det_corr.position, [], 3);
    detector.Nslice = size(det_corr.position, 2);
    detector.startslice = 1;
    detector.endslice = detector_corr.Nslice;
    detector.hx_ISO = det_corr.hx_ISO;
    detector.hz_ISO = det_corr.hz_ISO;
    detector.slicemerge = 1:detector.Nslice;
    detector.mergescale = 1;
else
    error(['Unknown collimator protocol ' collimator]);
end
end


function detector = hardcodeexposure(collimator, detector, det_corr)
% explain the collimator (hard code)

switch lower(collimator)
    case {'all', 'open'}
        % all on
        detector.position = reshape(det_corr.position, [], 3);
        detector.Nslice = size(det_corr.position, 2);
        detector.startslice = 1;
        detector.endslice = detector.Nslice;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:detector.Nslice;
        detector.mergescale = 1;
    case {'16x0.55', '16x0.5', 'u16 16x0.625'}
        sliceindex = 5:20;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 16;
        detector.startslice = 5;
        detector.endslice = 20;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:16;
        detector.mergescale = 1;
    case {'16x1.1', '16x1.0', 'u16 16x1.25'}
        detector.position = reshape(det_corr.position, [], 3);
        detector.Nslice = 24;
        detector.startslice = 1;
        detector.endslice = 24;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO*2;
        detector.slicemerge = [1 1 2 2 3 3 4 4 5 6 7 8 9 10 11 12 13 13 14 14 15 15 16 16];
        detector.mergescale = 2;
    case {'8x1.1', '8x1.0', 'u16 8x1.25'}
        sliceindex = 5:20;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 16;
        detector.startslice = 5;
        detector.endslice = 20;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO*2;
        detector.slicemerge = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
        detector.mergescale = 2;
    case {'8x0.55', '8x0.5', 'u16 8x0.625'}
        sliceindex = 8:15;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 8;
        detector.startslice = 8;
        detector.endslice = 15;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:8;
        detector.mergescale = 1;
    case {'4x0.55', '4x0.5', 'u16 4x0.625'}
        sliceindex = 10:13;
        detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
        detector.Nslice = 4;
        detector.startslice = 10;
        detector.endslice = 13;
        detector.hx_ISO = det_corr.hx_ISO;
        detector.hz_ISO = det_corr.hz_ISO;
        detector.slicemerge = 1:4;
        detector.mergescale = 1;
    otherwise
        error(['Unknown collimator protocol ', protocol.collimator]);  
end
end