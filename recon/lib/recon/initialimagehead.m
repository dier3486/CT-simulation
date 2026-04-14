function imagehead = initialimagehead(poolsize, imageheadfields)

if nargin<1 || ~isavail(poolsize)
    poolsize = 0;
end
if nargin<2
    % no image head?
    imageheadfields = {};
end

imagehead = struct();
if isincell(imageheadfields, 'InstanceNumber')
    % DICOM tag, index of the image, start from 1
    imagehead.InstanceNumber = zeros(1, poolsize, 'int32');
end
if isincell(imageheadfields, 'imagecenter')
    % XYZ position of the imagecenter to ISO
    imagehead.imagecenter = zeros(3, poolsize, 'single');
end
if isincell(imageheadfields, 'reconcenter')
    % reconcenter_2DBP
    imagehead.reconcenter = zeros(3, poolsize, 'single');
    % the center of the ring artifacts
end
if isincell(imageheadfields, 'ShotNumber')
    % shot index ?
    imagehead.ShotNumber = zeros(1, poolsize, 'int32');
end
if isincell(imageheadfields, 'SliceLocation')
    % DICOM tag, Z position of the imagecenter while couch on zero
    imagehead.SliceLocation = zeros(1, poolsize, 'single');
end
if isincell(imageheadfields, 'ImagePositionPatient')
    % DICOM tag, patient sight image position
    imagehead.ImagePositionPatient = zeros(3, poolsize, 'single');
end
if isincell(imageheadfields, 'Nnonzero')
    % nnz of the (sparse) image
    imagehead.Nnonzero = zeros(1, poolsize, 'int32');
end

% other head fields (should not have)
definedimagehead = {'InstanceNumber', 'imagecenter', 'reconcenter', 'ShotNumber', ...
        'SliceLocation', 'ImagePositionPatient', 'Nnonzero'};
unknowimagehead = setdiff(imageheadfields, definedimagehead);
if ~isempty(unknowimagehead)
    for ifield = unknowimagehead
        imagehead.(ifield{1}) = zeros(0, poolsize, 'single');
    end
end

end