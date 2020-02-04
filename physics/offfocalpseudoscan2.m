function Aoff = offfocalpseudoscan2(SYS, P, Pair)
% off-focal (pseudo) simulation base on convolution method (offfocalconv.m)
% only for Axial now
% Aoff = offfocalpseudoscan(SYS, P, Pair);
% e.g. Dataflow.P{iw} = Dataflow.P{iw}.*(1-offintensity)+offfocalpseudoscan2(SYS, Dataflow.P{iw}, Dataflow.Pair{iw});
% where the P{iw} is the output of projectionscan.m and reshaped like (Nps * Nview).

% only for Axial
if ~strcmpi(SYS.protocol.scan, 'Axial')
    warning('Off-focal simulation is not supported in %s scan now! sorry', SYS.protocol.scan);
    Aoff = [];
    return
end

% paramters to use
offwidth = SYS.source.offfocalwidth;
offintensity = SYS.source.offfocalintensity;
focalposition = SYS.source.focalposition;
Nfocal = SYS.source.focalnumber;
Npixel = double(SYS.detector.Npixel);
Nslice = double(SYS.detector.Nslice);
Nviewprot = SYS.protocol.viewperrot;
Nshot = SYS.protocol.shotnumber;
% if air
Nviewprot = min(Nviewprot, size(P, 2));

% slice weight
% w_slice = weightofslicemerge(SYS.detector);

% ini output
Aoff = zeros(size(P));

% fly focal
for ishot = 1:Nshot
    for ifocal = 1:Nfocal
        % get the views
        viewindex = (ifocal:Nfocal:Nviewprot) + (ishot-1)*Nviewprot;
        A_ii = P(:, viewindex)./Pair(:, ifocal);
        % mean on slice
        A_ii = mean(reshape(A_ii, Npixel, Nslice, []), 2);
        % off-focal convolution
        A_ii = offfocalconv(A_ii, SYS.detector, focalposition(ifocal, :), Nviewprot/Nfocal, offwidth, offintensity);
        % cross by w_slice and reshape
        A_ii = permute(reshape(repmat(A_ii(:), 1, Nslice), Npixel, Nviewprot/Nfocal, Nslice), [1 3 2]);
        A_ii = reshape(A_ii, Npixel*Nslice, []);
        % to return
        Aoff(:, viewindex) = A_ii.*Pair(:, ifocal);
    end
end

end