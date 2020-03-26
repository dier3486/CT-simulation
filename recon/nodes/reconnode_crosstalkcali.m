function [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)
% crosstalk calibration, odd-symetric style
% [dataflow, prmflow, status] = reconnode_crosstalkcali(dataflow, prmflow, status)

% parameters to use in prmflow
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;
Nps = Npixel*Nslice;
Nview = prmflow.recon.Nview;

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
datafields = findfields(dataflow, '\<rawdata_bk');
Nbk = length(datafields);
headfields = findfields(dataflow, '\<rawhead_bk');

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
    NLcorr.main = reshape(NLcorr.main, Nps, []);
else
    NLcorr = [];
end

% I know
HCscale = 1000;

% inverse the ideal data
for ibk = 2:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
    for ipx = 1:Nps
        % inverse non-linear corr
        if ~isempty(NLcorr)
            dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
        end
%         % inverse beamharden
%         dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
    end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end
% inverse the original data
for ibk = 1:2:Nbk
    % inverse Housefield
    dataflow.(datafields{ibk}) = dataflow.(datafields{ibk})./HCscale;
%     for ipx = 1:Nps
%         % apply the non-linear corr
%         dataflow.(datafields{ibk})(ipx, :) = iterpolyval(NLcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
%         % inverse beamharden
% %         dataflow.(datafields{ibk})(ipx, :) = iterinvpolyval(BHcorr.main(ipx, :), dataflow.(datafields{ibk})(ipx, :));
%     end
%     % inverse air (almost)
%     dataflow.(datafields{ibk}) = dataflow.(datafields{ibk}) + airrate;
end
% we should move above codes in another node specially do inverse

% % to intensity
% for ibk = 1:Nbk
%     dataflow.(datafields{ibk}) = 2.^(-dataflow.(datafields{ibk}));
% end

% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

p_order = 2;
alpha = 0.0002;
% Lcrs = [0 0 0; -1 2 -1; 0 0 0];
Lcrs = [1    -1     0
       -1     2    -1
        0    -1     1]./2;

L3 = [0 0 0; -1 2 -1; 0 0 0]./2;    
    
% loop slice
Pcrs = zeros(Npixel*Nslice, p_order);
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

%     Ap = spdiags([p(:,2) p(:,1)], [-1 1], Npixel,Npixel)';
    
    p = nan(Npixel, p_order);
    p2 = nan(Npixel, 2);
    p3 = nan(Npixel, 3);
    p4 = nan(Npixel, 3);
%     lamda = 0.0;
%     options = optimoptions('lsqnonlin','Display','off');
    % loop pixel
    for ipixel = 10:Npixel-9
        x = zeros(5, Nview*Nbk/2);
        y = zeros(5, Nview*Nbk/2);
        pixelindex = ipixel-2:ipixel+2;
        for ibk = 1:Nbk/2
            ix = ibk*2-1;
            iy = ibk*2;
            viewindex = (1:Nview) + (ibk-1)*Nview;
            x(:, viewindex) = double(dataflow.(datafields{ix})(pixelindex + Npixel*(islice-1), :));
            y(:, viewindex) = double(dataflow.(datafields{iy})(pixelindex + Npixel*(islice-1), :));
        end
        s = all(Seff(pixelindex, :), 1);
        if any(sum(reshape(s, Nbk/2, Nview),2)<Nview/2)
            continue;
        end
%         t0 = [0, 0];
%         p(ipixel, :) = lsqnonlin(@(t) crossfit3(t, x, y, s, rrate, lamda), t0, [], [], options);

        % test s
%         s(1152:end) = false;
%         s(1:1152) = false;
        w = 1./y(3, s);
%         w = 1.0;
        % p
        p(ipixel, :) = crossfit(x, y, s, w);
        % order 2
%         p2(ipixel, :) = crossfit(x, y, s, w);
%         % order 3
%         p3(ipixel, :) = crossfit3(x, y, s, w, L3.*alpha);
%         % order 4
%         p4(ipixel, :) = crossfit4(x, y, s, w);
        
%         if ipixel==416
%             1;
%         end
    end
    
%     % remove 0-filter
%     if p_order==3
%         d0 = mean(p(:,2), 'omitnan');
%         p = p + [d0/2 -d0 d0/2];
%     end
    
%     p = p - mean(sum(p,2), 'omitnan')/2;
    p = reshape(p, Npixelpermod, Nmod, p_order);
    s_unv = any(any(isnan(p),1),3);
    p(:, s_unv, :) = nan;
    pmod = mean(p, 2, 'omitnan');
    p(:, s_unv, :) = repmat(pmod, 1, sum(s_unv));
%     % try1, replace all module by mean
%     p = repmat(pmod, 1, Nmod);

    p = reshape(p, Npixel, p_order);
% %     % try2, smooth
% 	p(:,1) = p(:,1) - smooth(p(:,1), 0.04, 'loess');
%     p(:,2) = p(:,2) - smooth(p(:,2), 0.04, 'loess');
    
    p(1,1) = 0; p(end,end) = 0;
    
    pindex = (1:Npixel) + (islice-1)*Npixel;
    Pcrs(pindex, :) = p;
end

% slice merge
Pcrs = reshape(Pcrs, Npixel, Nmerge, []);
Pcrs = reshape(repmat(mean(Pcrs, 2), 1, Nmerge), [], p_order);

% % fill zero
if mod(p_order,2)==0
    Pcrs = [Pcrs(:,1:p_order/2) zeros(Nps, 1) Pcrs(:,p_order/2+1:end)];
end

% % m2, debug
% Pcross = zeros(Nps, 3);
% s = Pcrs(:,1)>0;
% Pcross(s, 2) = -Pcrs(s, 1).*2;
% Pcross(s, 1) = Pcrs(s, 1).*2;
% Pcross(~s, 2) = Pcrs(~s, 1).*2;
% Pcross(~s, 3) = -Pcrs(~s, 1).*2;

% paramters for corr
crosstalkcorr = caliprmforcorr(prmflow, corrversion);
% copy results to corr
crosstalkcorr.Nslice = Nslice;
crosstalkcorr.order = floor(p_order/2)*2+1;
crosstalkcorr.mainsize = Nps*crosstalkcorr.order;
crosstalkcorr.main = Pcrs;

% to return
dataflow.crosstalkcorr = crosstalkcorr;

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];

end


function p = crossfit(x, y, s, w)
% cross coefficients fitting function

m = size(x, 1);
mhf = (m+1)/2;

% w = 1./(-y(2,s).*log2(y(2,s)));
x13w = [x(mhf-1,s).*w; x(mhf+1,s).*w];

A = x13w * x13w';
b = x13w * ((y(mhf,s)-x(mhf,s)).*w)';
p = A\b;

p = p';

end

function p = crossfit3(x, y, s, w, alpha)
% cross coefficients fitting function

m = size(x, 1);
mhf = (m+1)/2;

xw = x(mhf-1:mhf+1, s).*w;
A = xw * xw';
A = A + diag(diag(A))*alpha;

b = xw * ((y(mhf,s)-x(mhf,s)).*w)';
p = A\b;

p = p';

end

function p = crossfit4(x, y, s, w)
% cross coefficients fitting function

m = size(x, 1);
mhf = (m+1)/2;
index = mhf + [-2 -1 1 2];

xw = x(index, s).*w;
A = xw * xw';

b = xw * ((y(mhf,s)-x(mhf,s)).*w)';
p = A\b;

p = p';

end