function Aoff = offfocalpseudoscan(SYS, P)
% off-focal (pseudo) simulation base on convolution method (offfocalconv.m)
% only for Axial now
% Aoff = offfocalpseudoscan(SYS, P);
% e.g. Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity)+offfocalpseudoscan(SYS, Dataflow.P{iw});
% where the P{iw} is the output of projectionscan.m and reshaped like (Nps * Nview).

% only for Axial
if ~strcmpi(SYS.protocol.scan, 'Axial')
    warning('Off-focal simulation is not supported in %s scan now! sorry', SYS.protocol.scan);
    Aoff = zeros(size(P));
    return
end

% paramters to use
offwidth = SYS.source.offfocalwidth;
offintensity = SYS.source.offfocalintensity;
if isfield(SYS.source, 'offedge')
    offedge = SYS.source.offedge;
else
    offedge = 1;
end
focalposition = SYS.source.focalposition;
Nfocal = SYS.source.focalnumber;
Npixel = double(SYS.detector.Npixel);
Nslice = double(SYS.detector.Nslice);
Nviewprot = SYS.protocol.viewperrot;
Nshot = SYS.protocol.shotnumber;
% if air
Nviewprot = min(Nviewprot, size(P, 2));

% slice weight
w_slice = weightofslicemerge(SYS.detector);

% ini output
Aoff = zeros(size(P));

% fly focal
for ishot = 1:Nshot
    for ifocal = 1:Nfocal
        % get the views
        viewindex = (ifocal:Nfocal:Nviewprot) + (ishot-1)*Nviewprot;
        A_ii = P(:, viewindex);
        % mean on slice
        A_ii = sum(reshape(A_ii, Npixel, Nslice, []), 2)./sum(w_slice);
        % off-focal convolution
        A_ii = offfocalconv(A_ii, SYS.detector, focalposition(ifocal, :), Nviewprot/Nfocal, offwidth, offintensity, offedge);
        % cross by w_slice and reshape
        A_ii = permute(reshape(A_ii(:)*w_slice, Npixel, Nviewprot/Nfocal, Nslice), [1 3 2]);
        A_ii = reshape(A_ii, Npixel*Nslice, []);
        % to return
        Aoff(:, viewindex) = A_ii;
    end
end

end