% pin cali main, test
% pinxml = 'F:\data-Dier.Z\PG\bay3_0601\Pinrecon.xml';
% pindata = 'F:\data-Dier.Z\PG\bay3_0601\pin\pin_1';
% pinxml = 'F:\data-Dier.Z\PG\bay4\Pinrecon.xml';
% pindata = 'F:\data-Dier.Z\PG\bay4\20200609\pin\2.1591612148003';
pinxml = 'F:\data-Dier.Z\PG\bay3\Pinrecon.xml';
pindata = 'F:\data-Dier.Z\PG\bay3\20200609\pin\smallfocal_airbowtie';
[~, dataflow, prmflow] = CRISrecon(pinxml, pindata);


% inputs are dataflow and prmflow
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nrot = Nview/Nviewprot;
detector = prmflow.system.detector;
focalposition = prmflow.system.focalposition(prmflow.recon.focalspot, :);
Npixelpermod = 16;
Nmod = Npixel/Npixelpermod;
SID = double(detector.SID);
hz = double(detector.hz_ISO);
viewangle = double(dataflow.rawhead.viewangle);

% det to fanangle
[fanangles0, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles0 = fanangles0 - focalangle;
fanangles0 = reshape(fanangles0, Npixel, Nslice);

% ini
    x0 = [200/SID  pi/2     0       0       0         0     0       0       0       0];
%        [r,       phi,     zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, viewoff, viewphanse,  zshift]
    s1 = [1        1        1       1       1         1     1       1       1       1];

alphaL = [1.5 0.1];
% fanfix
[fanfix, pinfit] = pindetcalibration(dataflow.rawdata, viewangle, detector, focalposition, Npixel, Nslice, Npixelpermod, Nviewprot, alphaL, x0);
% [fanfix, pinfit] = pindetcali_test1(dataflow, prmflow, x0);

midfanfix = pinfit(6);
focalpos = prmflow.system.focalposition(prmflow.recon.focalspot, :);
Xfix = (detector.position(:, 1) - focalpos(1)).*cos(fanfix(:) - midfanfix) - ...
       (detector.position(:, 2) - focalpos(2)).*sin(fanfix(:) - midfanfix) + focalpos(1);
Yfix = (detector.position(:, 1) - focalpos(1)).*sin(fanfix(:) - midfanfix) + ...
       (detector.position(:, 2) - focalpos(2)).*cos(fanfix(:) - midfanfix) + focalpos(2);
    
detector_fix = detector;
detector_fix.position(:, [1 2]) = [Xfix, Yfix];

% outputcorr = 'F:\data-Dier.Z\PG\bay4_0603\detector_PG_bay4_0603p2_v1.0.corr';
% corrcfg = readcfgfile(cfgmatchrule(outputcorr));
% packstruct(detector_fix, corrcfg, outputcorr);

