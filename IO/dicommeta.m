function dcminfo = dicommeta(protocol, clinicalinfo)
% default dicom info

if nargin>1 
    dcminfo = clinicalinfo;
else
    dcminfo = struct();
end

% SourceApplicationEntityTitle
% if known

% SpecificCharacterSet
dcminfo.SpecificCharacterSet = 'GB18030';

% ImageType
dcminfo.ImageType = 'ORIGINAL\PRIMARY\AXIAL';

% InstanceCreationDate
dcminfo.InstanceCreationDate = datestr(now,'yyyymmdd');

% InstanceCreationTime
dcminfo.InstanceCreationTime = datestr(now, 'hhMMss');

% StudyDate, StudyTime (and others)
% if known

% Modality
dcminfo.Modality = 'CT';

% Manufacturer
dcminfo.Manufacturer = 'CT-simulation';

% ConversionType
dcminfo.ConversionType = 'WSD';

% ReferringPhysicianName
% if konwn

% StudyDescription
if isfield(protocol, 'StudyDescription')
    dcminfo.StudyDescription = protocol.StudyDescription;
end

% SeriesDescription
if isfield(protocol, 'SeriesDescription')
    dcminfo.SeriesDescription = protocol.SeriesDescription;
end

% PatientName, PatientID and 
% if knowm

% BodyPartExamined
if isfield(protocol, 'BodyPartExamined')
    dcminfo.BodyPartExamined = protocol.BodyPartExamined;
end

% SliceThickness
if isfield(protocol, 'imagethickness')
    dcminfo.SliceThickness = protocol.imagethickness;
end
% image slice thichness

% KVP
if isfield(protocol, 'KV')
    dcminfo.KVP = protocol.KV;
end

% DataCollectionDiameter
% maxFOV
if isfield(protocol, 'maxFOV')
    dcminfo.DataCollectionDiameter = protocol.maxFOV;
end

% SoftwareVersions
% to be

% ProtocolName
if isfield(protocol, 'ProtocolName')
    dcminfo.ProtocolName = protocol.ProtocolName;
end

% ReconstructionDiameter
% reconFOV
if isfield(protocol, 'reconFOV')
    dcminfo.ReconstructionDiameter = protocol.reconFOV;
end

% DistanceSourceToPatient
% SID

% GantryDetectorTilt
if isfield(protocol, 'gantrytilt')
    dcminfo.GantryDetectorTilt = protocol.gantrytilt;
end

% TableHeight
% couchheight + 0

% RotationDirection
dcminfo.RotationDirection = 'CW';
% fixed on clockwise

% ExposureTime
% to be, for Axial it = protocol.collimation*protocol.rotationspeed/protocol.shotcouchstep

% XRayTubeCurrent
if isfield(protocol, 'mA')
    dcminfo.XRayTubeCurrent = protocol.mA;
end

% Exposure
% to be, mAs, for Axial it = protocol.mA*protocol.rotationspeed

% GeneratorPower
if isfield(protocol, 'mA') && isfield(protocol, 'KV')
    dcminfo.GeneratorPower = protocol.mA*protocol.KV;
end
% kW or W ??

% RevolutionTime
if isfield(protocol, 'RevolutionTime')
    dcminfo.RevolutionTime = protocol.RevolutionTime;
elseif isfield(protocol, 'rotationspeed')
    dcminfo.RevolutionTime = protocol.rotationspeed;
    % I know it was named by 'rotationspeed'
end

% FocalSpots
if isfield(protocol, 'focalwidth')
    dcminfo.FocalSpots = protocol.focalwidth(1);
end
% It is not the protocol.focalsize

% ConvolutionKernel
if isfield(protocol, 'reconkernel')
    dcminfo.ConvolutionKernel = protocol.reconkernel;
end

% PatientPosition
if isfield(protocol, 'PatientPosition')
    dcminfo.PatientPosition = protocol.PatientPosition;
elseif ~isfield(dcminfo, 'PatientPosition')
    dcminfo.PatientPosition = 'HFS';
end

% Tube Angle
if isfield(protocol, 'TubeAngle')
    dcminfo.TubeAngle = protocol.TubeAngle;
end

% SingleCollimationWidth
if isfield(protocol, 'SingleCollimationWidth')
    dcminfo.SingleCollimationWidth = protocol.SingleCollimationWidth;
elseif isfield(protocol, 'CollimatedSliceThickness')
    dcminfo.SingleCollimationWidth = protocol.CollimatedSliceThickness;
end
% which is the nominal slicethickness (of collimation)

% TotalCollimationWidth
if isfield(protocol, 'TotalCollimationWidth')
    dcminfo.TotalCollimationWidth = protocol.TotalCollimationWidth;
end

% TableSpeed
if isfield(protocol, 'couchspeed')
    dcminfo.TableSpeed = protocol.couchspeed;
elseif isfield(protocol, 'TableSpeed')
    dcminfo.TableSpeed = protocol.TableSpeed;
end

% Table Feed per Rotation
if isfield(protocol, 'TableFeedperRotation')
    dcminfo.TableFeedperRotation = protocol.TableFeedperRotation;
end

% SpiralPitchFactor
% Ratio of the Table Feed per Rotation (0018,9310) to the Total Collimation Width (0018,9307).
if isfield(protocol, 'SpiralPitchFactor')
    dcminfo.SpiralPitchFactor = protocol.SpiralPitchFactor;
end

% Data Collection Center (Patient)
% (0018,9313) The x, y, and z coordinates (in the Patient-Based Coordinate System) in mm of the center of the region in which
% data were collected. 
if isfield(protocol, 'DataCollectionCenterPatient')
    dcminfo.DataCollectionCenterPatient = protocol.DataCollectionCenterPatient;
end

% Reconstruction Target Center (Patient)
% (0018,9318) The x, y, and z coordinates (in the Patient-Based Coordinate System) of the reconstruction center target point 
% as used for reconstruction in mm.
if isfield(protocol, 'ReconstructionTargetCenterPatient')
    dcminfo.ReconstructionTargetCenterPatient = protocol.ReconstructionTargetCenterPatient;
end

% StudyInstanceUID
if isfield(protocol, 'StudyInstanceUID')
    dcminfo.StudyInstanceUID = protocol.StudyInstanceUID;
end

% SeriesInstanceUID
if isfield(protocol, 'SeriesInstanceUID')
    dcminfo.SeriesInstanceUID = protocol.SeriesInstanceUID;
end

% StudyID
if isfield(protocol, 'StudyID')
    dcminfo.StudyID = protocol.StudyID;
elseif ~isfield(dcminfo, 'StudyID')
    dcminfo.StudyID = '0';
end

% SeriesNumber
if isfield(protocol, 'seriesindex')
    dcminfo.SeriesNumber = protocol.seriesindex;
elseif isfield(protocol, 'SeriesNumber')
    dcminfo.SeriesNumber = protocol.SeriesNumber;
end

% AcquisitionNumber
if ~isfield(dcminfo, 'AcquisitionNumber')
    dcminfo.AcquisitionNumber = 1;
end

% Instance Number 
% shot index
% later

% Slice Thinkness
% The nominal slice thickness (of images)
if isfield(protocol, 'imagethickness')
    dcminfo.SliceThinkness = protocol.imagethickness;
end
% I know the protocol.imagethickness is nominal thickness and the 'real'
% thickness is recon.imagethickness

% Image Position (Patient)
% The x, y, and z coordinates of the upper left hand corner (center of the first voxel transmitted) of the image, in mm. 
if isfield(protocol, 'ImagePositionPatient')
    dcminfo.ImagePositionPatient = protocol.ImagePositionPatient;
end

% Slice Location
% later

% Pixel Spacing
if isfield(protocol, 'PixelSpacing')
    dcminfo.PixelSpacing = protocol.PixelSpacing;
end

% SamplesPerPixel
dcminfo.SamplesPerPixel = 1;

% PhotometricInterpretation
% default

% Rows & Columns
% AUTO

% BitsAllocated & BitsStored
% AUTO

% HighBit
% AUTO

% PixelRepresentation 
% AUTO

% SmallestImagePixelValue & LargestImagePixelValue
% AUTO

% Image Orientation (Patient)
% The direction cosines of the first row and the first column with respect to the patient.
if isfield(protocol, 'ImageOrientationPatient')
    dcminfo.ImageOrientationPatient = protocol.ImageOrientationPatient;
end
% I know for HFS tilt, it = [1 0 0  0 cos(tilt) sin(tilt)];

% FrameOfReferenceUID
if ~isfield(dcminfo, 'FrameOfReferenceUID')
    dcminfo.FrameOfReferenceUID = '';
end

% PositionReferenceIndicator
% no?
if ~isfield(dcminfo, 'PositionReferenceIndicator')
    dcminfo.PositionReferenceIndicator = '';
end

% WindowCenter
if isfield(protocol, 'windowcenter')
    dcminfo.WindowCenter = protocol.windowcenter;
end

% WindowWidth
if isfield(protocol, 'windowwidth')
    dcminfo.WindowWidth = protocol.windowwidth;
end

% RescaleIntercept & RescaleSlope & RescaleType
% fixed
dcminfo.RescaleIntercept = -1200;
dcminfo.RescaleSlope = 0.125;
dcminfo.RescaleType = 'HU';


end