function detector = detectorextracollim(detector, det_corr, samplekeV)
% merge extra informations of detector by collimation
% detector = detectorextracollim(detector, det_corr, samplekeV);

% response
if isfield(det_corr, 'response') && ~isempty(det_corr.response)
    % samplekeV of the response curve(s)
    detector.samplekeV = det_corr.samplekeV;
    % collimation
    if size(det_corr.response, 1)==1
        detector.response = det_corr.response;
    elseif size(det_corr.response, 1)==det_corr.Nslice
        % different slice different response (a simplified model of kneel effect)
        det_corr.response = repmat(reshape(det_corr.response, 1, det_corr.Nslice, []), det_corr.Npixel, 1, 1);
        sliceindex = detector.startslice : detector.endslice;
        detector.response = reshape(det_corr.response(:, sliceindex, :), detector.Npixel*detector.Nslice, []);
    else
        det_corr.response = reshape(det_corr.response, det_corr.Npixel, det_corr.Nslice, []);
        sliceindex = detector.startslice : detector.endslice;
        detector.response = reshape(det_corr.response(:, sliceindex, :), detector.Npixel*detector.Nslice, []);
    end
    % interp to samplekeV
    [index, alpha] = interpprepare(detector.samplekeV, samplekeV, 0);
    detector.response = detector.response(:, index(:,1)).*alpha(:,1)' + detector.response(:, index(:,2)).*alpha(:,2)';
else
    detector.response = ones(size(samplekeV));
end

% cross talk
if isfield(det_corr, 'crossmatrix') && ~isempty(det_corr.crossmatrix)
    % crossmatrix
    crsindex = false(det_corr.Npixel, det_corr.Nslice);
    crsindex(:, detector.startslice : detector.endslice) = true;
    detector.crossmatrix = det_corr.crossmatrix(crsindex(:), crsindex(:));
end

% norm vector
if isfield(det_corr, 'normvector') && ~isempty(det_corr.normvector)
    det_corr.normvector = reshape(det_corr.normvector, det_corr.Npixel, det_corr.Nslice, []);
    sliceindex = detector.startslice : detector.endslice;
    detector.normvector = reshape(det_corr.normvector(:, sliceindex, :), [], 3);
else
    % default vector
    focalposition = reshape(det_corr.focalposition, [], 3);
    detector.normvector = focalposition(1, :) - detector.position;
    detector.normvector(:, 3) = 0;
    detector.normvector = normr(detector.normvector);
end

% edgelength & pixel area
if isfield(det_corr, 'edgelength') && ~isempty(det_corr.edgelength)
    % defined edgelength
    if size(det_corr.edgelength(:), 1) <= 2
        det_corr.edgelength = repmat(reshape(det_corr.edgelength, 1,1,[]), det_corr.Npixel, det_corr.Nslice, 1);
    else
        det_corr.edgelength = reshape(det_corr.edgelength, det_corr.Npixel, det_corr.Nslice, []);
    end
    sliceindex = detector.startslice : detector.endslice;
    detector.edgelength = det_corr.edgelength(:, sliceindex, :);
    % pixelarea = edgelength(1)*edgelength(2)
    detector.pixelarea = detector.edgelength(:, :, 1).*detector.edgelength(:, :, 2);
elseif isfield(det_corr, 'pixelarea') && ~isempty(det_corr.pixelarea)
    % defined pixelarea
    if size(det_corr.pixelarea(:), 1) == 1
        det_corr.pixelarea = ones(det_corr.Npixel, det_corr.Nslice).*det_corr.pixelarea;
    else
        det_corr.pixelarea = reshape(det_corr.pixelarea, det_corr.Npixel, det_corr.Nslice);
    end
    sliceindex = detector.startslice : detector.endslice;
    detector.pixelarea = det_corr.pixelarea(:, sliceindex);
    % default detector.edgelength(1)=1
    detector.edgelength = ones(detector.Npixel, detector.Nslice, 2);
    detector.edgelength(:,:,2) = detector.pixelarea./detector.edgelength(:,:,1);
else
    % default pixel area
    w = weightofslicemerge(detector);
    detector.pixelarea = repmat(w(:)', detector.Npixel, 1);
    % default detector.edgelength(1)=1
    detector.edgelength = ones(detector.Npixel, detector.Nslice, 2);
    detector.edgelength(:,:,2) = detector.pixelarea./detector.edgelength(:,:,1);
end
detector.pixelarea = detector.pixelarea(:);
detector.edgelength = reshape(detector.edgelength, [], 2);

% TBC
% ASG,

end
