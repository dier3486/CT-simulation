function detector = collimatorexposure(collimator, detector, det_corr, collimatorexplain)
% explain the collimator
% detector = collimatorexposure(collimator, detector, det_corr, collimatorexplain);

% prepare
det_corr.position = reshape(det_corr.position, det_corr.Npixel, det_corr.Nslice, 3);
detector.Npixel = det_corr.Npixel;
detector.SID = det_corr.SID;
detector.SDD = det_corr.SDD;

% no collimatorexplain?
if nargin<4
    collimatorexplain = [];
end

if isempty(collimatorexplain)
    % load default explain (hard code)
    detector = hardcodeexposure(collimator, detector, det_corr);
else
    % explain the collimator
    Nc = length(collimatorexplain(:));
    for icoll = 1:Nc
        coll = collimatorexplain{icoll};
        if strcmpi(collimator, coll.name)
            % get the collimator protocol
            Nslice_orig = size(det_corr.position, 2);
            detector.Nslice = coll.Nslice;
            % startslice and endslice
            if isfield(coll, 'startslice') && isfield(coll, 'endslice')
                detector.startslice = coll.startslice;
                detector.endslice = coll.endslice;
            else
                detector.startslice = (Nslice_orig - coll.Nslice)/2 + 1;
                detector.endslice = (Nslice_orig + coll.Nslice)/2;
            end
            % detector position
            sliceindex = detector.startslice : detector.endslice;
            detector.position = reshape(det_corr.position(:, sliceindex, :), [], 3);
            % slicemerge and/or mergescale
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
            % h
            detector.hx_ISO = det_corr.hx_ISO;
            detector.hz_ISO = det_corr.hz_ISO * detector.mergescale;
            
            % rowcombine
            if isfield(coll, 'rowcombine')
                detector.rowcombine = coll.rowcombine;
                % NOTE: empty coll.rowcombine means do nothing in recon node rowcombine.
            else
                % default rowcombine (hard code)
                detector.rowcombine = hardcoderowcombine(coll);
                % NOTE: you should set them in configure files.
            end
            % break
            return
        end
    end

    % why go to here?
    error(['Unknown collimator protocol ' collimator]);
%     if strcmpi(collimator, 'all') || strcmpi(collimator, 'open')
%         % all on for lazy
%         detector.position = reshape(det_corr.position, [], 3);
%         detector.Nslice = size(det_corr.position, 2);
%         detector.startslice = 1;
%         detector.endslice = detector_corr.Nslice;
%         detector.hx_ISO = det_corr.hx_ISO;
%         detector.hz_ISO = det_corr.hz_ISO;
%         detector.slicemerge = 1:detector.Nslice;
%         detector.mergescale = 1;
%         detector.rowcombine = [];
%     else
%     error(['Unknown collimator protocol ' collimator]);
%     end
end

% returned previously, do NOT put anything here
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
        detector.slicemerge = [1 2 3 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 14 15 16];
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
        error(['Unknown collimator protocol ', collimator]);  
end
detector.rowcombine = [];

end


function rowcombine = hardcoderowcombine(coll)
% default rowcombine (hardcode), only for 

colltoken = regexp(coll.name, '\<(\d+)[x \*]', 'tokens');
if isempty(colltoken)
    % what??
    rowcombine = [];
    return;
end
rownominal = str2double(colltoken{1}{1});

rowscale = coll.Nslice/rownominal;
if rowscale==floor(rowscale)
    rowcombine = repelem(1:rownominal, 1, rowscale);
else
    % hard code for PX
    switch lower(coll.name)
        case {'16x1.2'}
            rowcombine = [1 2 3 4 5 5 6 6 7 7 8 8 9 9 10 10 11 11 12 12 13 14 15 16];
        case {'8x1.2'}
            rowcombine = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8];
        otherwise
            % unknown collimator
            rowcombine = [];
    end
end
end
