function [fanfix, pinfit] = pindetcali_test1(dataflow, prmflow, x0)

% inputs are dataflow and prmflow
Nslice = prmflow.recon.Nslice;
Npixel = prmflow.recon.Npixel;
Nview = prmflow.recon.Nview;
Nviewprot = prmflow.recon.Nviewprot;
Nrot = Nview/Nviewprot;
detector = prmflow.system.detector;
focalposition = prmflow.system.focalposition(prmflow.recon.focalspot, :);
Npixelpermod = 16;
% Nmod = Npixel/Npixelpermod;
SID = double(detector.SID);
SDD = double(detector.SDD);
hz = double(detector.hz_ISO);
viewangle = double(dataflow.rawhead.viewangle);

if nargin<3
    x0 = [200/SID  0        0       0       0         0     0       0       0       0];
%        [r,       phi,     zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, viewoff, viewphanse,  zshift]
end

% det to fanangle
[fanangles0, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles0 = fanangles0 - focalangle;
fanangles0 = reshape(fanangles0, Npixel, Nslice);

% ini
fanangles = fanangles0;

figid = figure;
% subplot(2,1,1);
subplot(2,1,2);
hold on
% iteration of det fix
iter_stop = [2.0e-3, 2.0e-4];
alpha = 1.0;
Nitermax = 20;
pinfit = x0;
detfix0 = zeros(Npixel*Nslice, Nitermax+1);
for iiter = 1:Nitermax
    % fit splines for each view's projection
    cs = pinsplines(dataflow.rawdata, fanangles, Npixel, Nslice, Nview);
    % pin matrix
    alpha_L = [1.0, 0.1];
    [A, L, indexrange] = pinleftoptmatrix(cs, Npixel, Nslice, Nview, Npixelpermod, alpha_L);
    
    % pin curve fitting
%     p0 = reshape(double([cs(:).p]), Nslice, Nview);
    
    % fit pin
    [dp, pinfit] = pincurvefit(cs, viewangle, hz, Nslice, Nview, pinfit);
    
    detfix = zeros(Npixel, Nslice);
    d_right = zeros(Npixel, Nslice);
    for islice = 1:Nslice
        AA = A{islice}'*A{islice};
        edge1 = indexrange(islice, 1);  edge2 = indexrange(islice, 2);
        AA = AA(edge1:edge2, edge1:edge2) + L(edge1:edge2, edge1:edge2).*(Nrot^1.5);
        d_right(:, islice) = A{islice}'*dp(islice,:)';
        detfix(edge1:edge2, islice) = AA\d_right(edge1:edge2, islice);
    end
    % rec
    detfix0(:, iiter+1) = detfix(:);
    % normr
    norm_detfix = sqrt(mean(detfix(:).^2)) * SDD;
    norm_ddf = sqrt(mean((detfix(:) - detfix0(:, iiter)).^2)) * SDD;
    norminf_detfix = max(abs(detfix(:))) * SDD;
    norminf_ddf = max(abs(detfix(:) - detfix0(:, iiter))) * SDD;
    
    % print
    fprintf('%d: %e   %e   %e   %e\n', iiter, norm_detfix, norm_ddf, norminf_detfix, norminf_ddf);
    % fix det
    fanangles = fanangles + detfix.*alpha;
    % will be return
    fanfix = fanangles - fanangles0 - detfix.*alpha.*iiter;
    % plot
    if iiter>1
        figure(figid);
        midfix = pinfit(6)*SDD;
        subplot(2,1,1);  imagesc(fanfix'.*SDD - midfix, [-0.15 0.15] - midfix); colormap jet; colorbar;
        subplot(2,1,2);  plot(mean(detfix, 2));
        drawnow;
    end
    
    % is done?
    if norm_detfix<iter_stop(1) || norm_ddf<iter_stop(2)
        % done
        break;
    end
end

end


function [dp, pinfit] = pincurvefit(cs, viewangle, hz, Nslice, Nview, x0, s)

if nargin<7
    s = true(size(x0));
end
p0 = reshape(double([cs(:).p]), Nslice, Nview);
pinfit = lsqnonlin(@(x) pinfitfun(viewangle, Nslice, hz, x, p0, s), x0);
p1 = pinfitfun(viewangle, Nslice, hz, pinfit, 0, s);
p1 = p1(1:Nslice, 1:Nview);
dp = p1 - p0;

end


function r = pinfitfun(viewangle, Nslice, hz, x, p0, s)

x(s) = x;
x = num2cell(x);

p = pinprojectfun(viewangle, Nslice, hz, x{:});

r = p - p0;
er = (r - mean(r)).*(Nslice-1);
dr = r(:,end) - r(:, 1);
r = [r er dr.*size(p, 2)];

end