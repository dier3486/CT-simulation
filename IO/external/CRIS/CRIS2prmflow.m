function prmflow = CRIS2prmflow(prmflow, rawxml)
% put the parameters from rawxml to prmflow
% prmflow = CRIS2prmflow(prmflow, rawxml)

protocol_raw = struct();
% scan
protocol_raw.scan = rawxml.AcquisitionParameter.ScanType;
% collimator
protocol_raw.collimator = [num2str(rawxml.AcquisitionParameter.CollimatedSliceNum) 'x' ...
    num2str(rawxml.AcquisitionParameter.CollimatedSliceThickness)];
% bowtie
protocol_raw.bowtie = rawxml.AcquisitionParameter.bowTie;
% switch lower(rawxml.AcquisitionParameter.bowTie)
%     case 'air'
%         protocol_raw.bowtie = 'Empty';
%     case 'small'
%         protocol_raw.bowtie = 'Head';
%     case 'large'
%         protocol_raw.bowtie = 'Body';
%     otherwise
%         protocol_raw.bowtie = rawxml.AcquisitionParameter.bowTie;
% end
% focalspot
protocol_raw.focalspot = rawxml.AcquisitionParameter.FocusSpotSequence;
% focalsize
protocol_raw.focalsize = rawxml.AcquisitionParameter.FocusSpotSize;
% KV
protocol_raw.KV = rawxml.AcquisitionParameter.kVp;
% mA
protocol_raw.mA = rawxml.AcquisitionParameter.mA;
% view per rotation
protocol_raw.viewperrot = rawxml.AcquisitionParameter.ViewPerRevolution;
% rotationspeed
protocol_raw.rotationspeed = rawxml.AcquisitionParameter.gantrySpeed;
% rotationnumber (no)
% protocol_raw.rotationnumber = nan;
% viewnumber
protocol_raw.viewnumber = rawxml.AcquisitionParameter.numView;
% startangle
protocol_raw.startangle = rawxml.ReconParameters.startAngle;
% startcouch ?
protocol_raw.startcouch = rawxml.ReconParameters.startLoc;
% shotnumber
protocol_raw.shotnumber = rawxml.ReconParameters.EndShotNum - rawxml.ReconParameters.StartShotNum + 1;
% startshot
protocol_raw.startshot = rawxml.ReconParameters.StartShotNum;
% shotcouchstep (no?)
protocol_raw.shotcouchstep = rawxml.AcquisitionParameter.CollimatedSliceNum * ...
    rawxml.AcquisitionParameter.CollimatedSliceThickness;
% couchheight (no)
% protocol_raw.couchheight = 0;
% couchspeed
protocol_raw.couchspeed = rawxml.AcquisitionParameter.tableSpeed;
% pitch
protocol_raw.pitchhelical = rawxml.AcquisitionParameter.pitchHelical;
% recon kernel
protocol_raw.reconkernel = rawxml.ReconParameters.reconKernel;
% reconFOV
protocol_raw.reconFOV = rawxml.ReconParameters.displayFOV;
% ImageThickness
protocol_raw.imagethickness = rawxml.ReconParameters.ImageThickness;
% ImageIncrement
protocol_raw.imageincrement = rawxml.ReconParameters.ImageIncrement;
% imagesize
protocol_raw.imagesize = rawxml.ReconParameters.imgColumn;

% merge
if isfield(prmflow, 'protocol')
    prmflow.protocol = structmerge(prmflow.protocol, protocol_raw);
else
    prmflow.protocol = protocol_raw;
end

end