function [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)
% crosstalk calibration, odd-symetric style
% [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;
Nfocal = prmflow.recon.Nfocal;

% parameters to use
caliprm = prmflow.pipe.(status.nodename);
% caliprm = [];

% pixel number per module
if isfield(caliprm, 'Npixelpermod')
    Npixelpermod = caliprm.Npixelpermod;
elseif isfield(prmflow.system.detector, 'Npixelpermod')
    Npixelpermod = prmflow.system.detector.Npixelpermod;
else
    % default, 16 pixels in one module on x-direction
    Npixelpermod = 16;
end
Nmod = Npixel/Npixelpermod;
% slice merge in expanding
if isfield(caliprm, 'Nmerge')
    Nmerge = caliprm.Nmerge;
else
    % defualt to merge 4 slices
    Nmerge = min(4, Nslice);
end
% cross talk in intensity or log2
if isfield(caliprm, 'istointensity')
    istointensity = caliprm.istointensity;
else
    % defualt is to crosstalk in intensity
    istointensity = true;
end

% format version of calibration table
if isfield(caliprm, 'corrversion')
    corrversion = caliprm.corrversion;
else
    corrversion = 'v1.0';
end
% echo .
echo_onoff = true;

% inverse bh and nl
% find out the input data, I know they are rawdata_bk1, rawdata_bk2 ...
datafields = findfields(dataflow, '^rawdata_bk');
Nbk = length(datafields);
headfields = findfields(dataflow, '^rawhead_bk');

% in which odd are the original value, even are the ideal value (fitting target)
% reshape
for ibk = 1:Nbk
    dataflow.(datafields{ibk}) = reshape(dataflow.(datafields{ibk}), Nps, Nview);
end

% I know the air corr is in prmflow.corrtable
Aircorr = prmflow.corrtable.Air;
Aircorr.main = reshape(Aircorr.main, Nps, []);
% I know the beamharden corr is in prmflow.corrtable
BHcorr = prmflow.corrtable.Beamharden;
BHcorr.main = reshape(BHcorr.main, Nps, []);
% I know the nonlinear corr is in dataflow
if isfield(dataflow, 'nonlinearcorr')
    NLcorr = dataflow.nonlinearcorr;
    NLcorr.main = reshape(NLcorr.main, Nps, [], NLcorr.focalnumber);
    % NLcorr.focalnumber shall <= Nfocal
elseif isfield(prmflow.corrtable, 'Nonlinear')
    NLcorr = prmflow.corrtable.Nonlinear;
    NLcorr.main = reshape(NLcorr.main, Nps, [], NLcorr.focalnumber);
else
    NLcorr = [];
end

% I know
HCscale = 1000;

% inverse the ideal data
for ibk = 2:2:Nbk
    % inverse Hounsefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % inverse non-linear corr
        if ~isempty(NLcorr)
            for ifocal = 1:NLcorr.focalnumber
                dataflow.(datafields{ibk})(ipx, ifocal:NLcorr.focalnumber:end) = ...
                    iterinvpolyval( NLcorr.main(ipx, :, ifocal), ...
                    dataflow.(datafields{ibk})(ipx, ifocal:NLcorr.focalnumber:end));
            end
        end
        % inverse beamharden
        if istointensity
            dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        end
    end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end
% inverse the original data
for ibk = 1:2:Nbk
    % inverse Hounsefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % apply the non-linear corr
%         dataflow.(datafields{ibk})(ipx, :) = iterpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        % inverse beamharden
        if istointensity
            dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        end
    end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end
% we should move above codes in another node specially do inverse

% to intensity
if istointensity
    for ibk = 1:Nbk
        dataflow.(datafields{ibk}) = 2.^(-dataflow.(datafields{ibk}));
    end
end

% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

p_order = 1;    % only order 1 supported yet
alpha = 2.0;

% loop slice
Pcrs = zeros(Npixel*Nslice, p_order, Nfocal);
for islice = 1:Nslice
    if echo_onoff, fprintf('.'); end
    % Seff
    Seff = false(Npixel, Nview*Nbk/2);
    for ibk = 1:Nbk/2
        for iview=1:Nview
            viewindex = iview+(ibk-1)*Nview;
            index_ii = index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk);
            Seff(index_ii, viewindex) = true;
        end
    end
    pixavail = sum(Seff, 2) > (Nview*Nbk/2 /3);
    
    % get the samples
    x = zeros(Npixel, Nview*Nbk/2);
    y = zeros(Npixel, Nview*Nbk/2);
    for ibk = 1:Nbk/2
        ix = ibk*2-1;
        iy = ibk*2;
        viewindex = (1:Nview) + (ibk-1)*Nview;
        x(:, viewindex) = double(dataflow.(datafields{ix})((1:Npixel) + Npixel*(islice-1), :));
        y(:, viewindex) = double(dataflow.(datafields{iy})((1:Npixel) + Npixel*(islice-1), :));
    end
    % diff
    dx = [zeros(1, Nview*Nbk/2); diff(x);];
    dy = [zeros(1, Nview*Nbk/2); diff(y);];
    dx(~Seff) = 0;
    dy(~Seff) = 0;
    
    % The correction is fix the x_{i,j} by new x"_{i,j} = x_{i,j} + p_i*(x_{i+1,j}-x_{i-1,j}), where i is the pixel index 
    % and j is the view index, to make the x"_{i,j}-(x"_{i+1,j}+x"_{i-1,j})/2 approaching to y_{i,j}-(y_{i+1,j}+y_{i-1,j})/2.
    % Base on that the p_i is the solution of the these equations:
    
    % loop the focals
    p = nan(Npixel, p_order, Nfocal);
    for ifocal = 1:Nfocal 
        A = sparse([], [], [], Npixel, Npixel);
        b = zeros(Npixel, 1);
        for iview = ifocal: Nfocal: Nview*Nbk/2
            Aii = spdiags(([dx(2:end, iview); 0] + dx(1:end, iview))*[-1/2 1 -1/2], [-1 0 1], Npixel, Npixel);
%             Wii = spdiags(-log(y(:, iview)), 0, Npixel, Npixel);
%             Wii = spdiags(1./(-y(:,iview)./log(y(:,iview))), 0, Npixel, Npixel);
%             Wii = spdiags(1./y(:,iview), 0, Npixel, Npixel);
%             if iview<=Nview
%                 Wii = Wii.*2;
%             end
%             Wii = spdiags(1./y(:,iview).^2, 0, Npixel, Npixel);
            Wii = 1.0;
            A = A + Aii'*Wii*Aii;
            bii = (dy(:, iview) - [dy(2:end, iview); 0])./2 - (dx(:, iview) - [dx(2:end, iview); 0])./2;
            b = b + Aii'*Wii*bii;
        end
        b(~pixavail) = 0;
        L = speye(Npixel);
%         L = spdiags(ones(Npixel, 1)*[-1/2 1 -1/2], [-1 0 1], Npixel, Npixel);
        Lalpha = sqrt(sum(b.^2)./sum(pixavail)) * (alpha/0.1);
        % I konw crosstalk coefficients almost less than 0.1
        p(:, :, ifocal) = (A + L.*Lalpha)\b;
    end
    
    % remove unavailable pixels
    p(~pixavail) = nan;
    p = reshape(p, Npixelpermod, Nmod, p_order*Nfocal);
    s_unv = any(any(isnan(p),1),3);
    p(:, s_unv, :) = nan;
    pmod = mean(p, 2, 'omitnan');
    p(:, s_unv, :) = repmat(pmod, 1, sum(s_unv));
%     % try1, replace all module by mean
%     p = repmat(pmod, 1, Nmod);

    p = reshape(p, Npixel, p_order, Nfocal);
% %     % try2, smooth
% 	p(:,1) = p(:,1) - smooth(p(:,1), 0.04, 'loess');
%     p(:,2) = p(:,2) - smooth(p(:,2), 0.04, 'loess');
    
    p(1, 1, :) = 0; p(end, end, :) = 0;
    
    pindex = (1:Npixel) + (islice-1)*Npixel;
    Pcrs(pindex, :, :) = p;
end

% slice merge
Pcrs = reshape(Pcrs, Npixel, Nmerge, []);
Pcrs = reshape(repmat(mean(Pcrs, 2), 1, Nmerge), [], p_order*Nfocal);

% % fill zero
if mod(p_order,2)==0
%     Pcrs = [Pcrs(:,1:p_order/2) zeros(Nps, 1) Pcrs(:,p_order/2+1:end)];
    index_zz = true(1, (p_order+1)*Nfocal);
    index_zz(p_order/2+1 : p_order+1 : end) = false;
    Pcrs(:, index_zz) = Pcrs;
    Pcrs(:, ~index_zz) = 0;
end

% paramters for corr
crosstalkcorr = caliprmforcorr(prmflow, corrversion);
% copy results to corr
crosstalkcorr.Nslice = Nslice;
crosstalkcorr.order = floor(p_order/2)*2+1;
crosstalkcorr.istointensity = istointensity;
crosstalkcorr.mainsize = Nps*crosstalkcorr.order*Nfocal;
crosstalkcorr.main = Pcrs;

% to return
dataflow.crosstalkcorr = crosstalkcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end