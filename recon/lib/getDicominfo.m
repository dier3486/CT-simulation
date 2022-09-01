function dcminfo = getDicominfo(prmflow, status)
% get dicom info

% dicom meta
clinicalinfo = struct();
if isfield(prmflow, 'external')
    switch prmflow.external.Manufacturer
        case 'SINOVISION'
            % CRIS platform of SINOVISION
            clinicalinfo = CRISdicommeta(prmflow.protocol, prmflow.external.rawxml);
        otherwise
            1;
    end

end
dcminfo = dicommeta(prmflow.protocol, clinicalinfo);

% % debug
% dcminfo.ImageOrientationPatient = [1 0 0 0 cos(0.05) sin(0.05)]';
% dcminfo.GantryDetectorTilt = 0;

% UID
if ~isfield(dcminfo, 'StudyInstanceUID')
%     dcminfo.StudyInstanceUID = [status.taskUID '.' dcminfo.StudyID];
    dcminfo.StudyInstanceUID = status.taskUID;
end

if ~isfield(dcminfo, 'SeriesInstanceUID')
    dcminfo.SeriesInstanceUID = dicomuid();
end


% rep Nimage times
if isfield(prmflow.recon, 'Nimage')
    Nimage = prmflow.recon.Nimage;
else
    Nimage = 1;
end
dcminfo = repmat(dcminfo, Nimage, 1);

% set ImagePositionPatient and others
if isfield(prmflow.recon, 'imagecenter')
    XY = -prmflow.recon.FOV*(prmflow.recon.imagesize-1)/prmflow.recon.imagesize/2;
    for ii = 1:Nimage
        dcminfo(ii).ImagePositionPatient = prmflow.recon.imagecenter(ii, :) + [XY XY 0];
        dcminfo(ii).SliceLocation = prmflow.recon.imagecenter(ii, 3) + prmflow.recon.startcouch;
        dcminfo(ii).InstanceNumber = ii;
    end
else
    for ii = 1:Nimage
        dcminfo(ii).ImagePositionPatient = [0 0 ii];
        dcminfo(ii).SliceLocation = ii;
        dcminfo(ii).InstanceNumber = ii;
    end
end


end