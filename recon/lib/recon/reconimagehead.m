function imagehead = reconimagehead(imagehead, recon, imageindex, indexOut, ishot, Nnonzero)
% get imagehead

if nargin<1 || isempty(imagehead)
    % initial imagehead
    imagehead = initialimagehead();
end
if nargin<2
    return;
end
if nargin<4 || isempty(indexOut)
    indexOut = imageindex;
end
if nargin<5
    ishot = 1;
end
if nargin<6
    Nnonzero = [];
end


% InstanceNumber (image index)
imagehead.InstanceNumber(indexOut) = imageindex + recon.InstanceStart;
% imagecenter
imagehead.imagecenter(:, indexOut) = recon.imagecenter(:, imageindex);
% reconcenter
if isfield(recon, 'reconcenter_2DBP')
    imagehead.reconcenter(:, indexOut) = recon.reconcenter_2DBP(:, imageindex, :);
else
    imagehead.reconcenter(:, indexOut) = recon.imagecenter(:, imageindex);
end

% ShotNumber
imagehead.ShotNumber(indexOut) = ishot;
% SliceLocation
imagehead.SliceLocation(indexOut) = recon.imagecenter(3, imageindex) + recon.startcouch;
% ImagePositionPatient
voxelsize = recon.voxelsize;
XY = -voxelsize.*(recon.imagesize-1)./2;
imagehead.ImagePositionPatient(:, indexOut) = recon.imagecenter(:, imageindex) + [XY(:); 0];
% The InstanceNumber, SliceLocation and ImagePositionPatient are public DICOM tags.

% Nnonzero
if ~isempty(Nnonzero)
    imagehead.Nnonzero(indexOut) = int32(Nnonzero);
end

end
