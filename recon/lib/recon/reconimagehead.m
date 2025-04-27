function imagehead = reconimagehead(imagehead, recon, imageindex, index_out, ishot)
% get imagehead

if nargin<1 || isempty(imagehead)
    % initial imagehead
    imagehead = initialimagehead();
end
if nargin<2
    return;
end
if nargin<4 || isempty(index_out)
    index_out = imageindex;
end
if nargin<5
    ishot = 1;
end

% InstanceNumber (image index)
imagehead.InstanceNumber(index_out) = imageindex + recon.InstanceStart;
% imagecenter
imagehead.imagecenter(:, index_out) = recon.imagecenter(imageindex, :)';
% reconcenter
if isfield(recon, 'reconcenter_2DBP')
    imagehead.reconcenter(:, index_out) = recon.reconcenter_2DBP(imageindex, :)';
else
    imagehead.reconcenter(:, index_out) = recon.imagecenter(imageindex, :)';
end
% Shot_Number
imagehead.Shot_Number(index_out) = ishot;
% SliceLocation
imagehead.SliceLocation(index_out) = recon.imagecenter(imageindex, 3)' + recon.startcouch;
% ImagePositionPatient
voxelsize = recon.voxelsize;
XY = -voxelsize.*(recon.imagesize-1)./2;
imagehead.ImagePositionPatient(:, index_out) = recon.imagecenter(imageindex, :)' + [XY 0]';

% The InstanceNumber, SliceLocation and ImagePositionPatient are public DICOM tags.

end
