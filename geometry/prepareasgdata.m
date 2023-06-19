function ASGdata = prepareasgdata(detector, samplekeV)
% prepare ASG data

if ~isfield(detector, 'ASGheight') || ~any(detector.ASGheight>0)
    ASGdata = [];
    return;
end

% ASGdata.nVasg = detector.nVasg;
ASGdata.nVedgex = detector.nVedgex;
ASGdata.nVedgez = detector.nVedgez;

ASGdata.nASGplx = normr(cross(detector.nVasg, detector.nVedgez));
ASGdata.nASGplz = normr(cross(detector.nVasg, detector.nVedgex));

ASGdata.normvector = detector.normvector;
ASGdata.ASGheight = detector.ASGheight(:)';
ASGdata.edgelength = detector.edgelength;
ASGdata.edgegap = detector.edgegap;
ASGdata.ASGthickness = detector.ASGthickness(:)';
ASGdata.gap_p = [(detector.edgegap(:, 1) + detector.edgegap(:, 2))./2  (detector.edgegap(:, 3) + detector.edgegap(:, 4))./2];
ASGdata.gap_n = [(detector.edgegap(:, 1) - detector.edgegap(:, 2))./2  (detector.edgegap(:, 3) - detector.edgegap(:, 4))./2];
ASGdata.pixelspace = detector.edgelength + ASGdata.gap_p.*2;

nAdotndet = dot(detector.nVasg, detector.normvector, 2);
ASGdata.nAdotnXondet = dot(detector.nVasg, ASGdata.nVedgex, 2)./nAdotndet;
ASGdata.nAdotnZondet = dot(detector.nVasg, ASGdata.nVedgez, 2)./nAdotndet;

% material
ASGdata.mu = interp1(detector.ASG.material.samplekeV, detector.ASG.material.mu_total, samplekeV);

end