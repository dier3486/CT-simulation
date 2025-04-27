function imagehead = initialimagehead(poolsize)

if nargin<1 || ~isavail(poolsize)
    poolsize = 0;
end

imagehead = struct();
if ~isfield(imagehead, 'InstanceNumber')
    % DICOM tag, index of the image, start from 1
    imagehead.InstanceNumber = zeros(1, poolsize, 'int32');
end
if ~isfield(imagehead, 'imagecenter')
    % XYZ position of the imagecenter to ISO
    imagehead.imagecenter = zeros(3, poolsize, 'single');
end
if ~isfield(imagehead, 'reconcenter')
    % reconcenter_2DBP
    imagehead.reconcenter = zeros(3, poolsize, 'single');
    % the center of the ring artifacts
end
if ~isfield(imagehead, 'Shot_Number')
    % shot index ?
    imagehead.Shot_Number = zeros(1, poolsize, 'int32');
end
if ~isfield(imagehead, 'SliceLocation')
    % DICOM tag, Z position of the imagecenter while couch on zero
    imagehead.SliceLocation = zeros(1, poolsize, 'single');
end
if ~isfield(imagehead, 'ImagePositionPatient')
    % DICOM tag, patient sight image position
    imagehead.ImagePositionPatient = zeros(3, poolsize, 'single');
end


end