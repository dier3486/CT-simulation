function [fanfix, pinfit] = pindetcalibration(rawdata, viewangle, detector, focalposition, Npixel, Nslice, Npixelpermod, ...
                            Nviewprot, alphaL, pinfit0, coeftofit)
% pin based detector postion calibration

Nview = size(viewangle(:), 1);
Nrot = Nview/Nviewprot;
% focalposition = prmflow.system.focalposition(prmflow.recon.focalspot, :);
SID = double(detector.SID);
SDD = double(detector.SDD);
hz = double(detector.hz_ISO);

if nargin<9
    alphaL = [1.0, 0.0];
end
if nargin<10
    pinfit0 = [200/SID  pi/2        0       0       0         0     0       0       0       0         0        0];
%             [r,       phi,     zeta_x, zeta_y, rotscale, midc, dslope, dscale, isooff, isophase, isooff2, isophase2,  zshift]
end
if nargin<11
    coeftofit = true(size(pinfit0));
else
    coeftofit = logical(coeftofit);
end

% det to fanangle
[fanangles0, focalangle] = detpos2fanangles(detector.position, focalposition);
fanangles0 = fanangles0 - focalangle;
fanangles0 = reshape(fanangles0, Npixel, Nslice);

% viewangle to double for matlab lsqnonlin
viewangle = double(viewangle);

% ini
fanangles = fanangles0;

figid = figure;
% subplot(2,1,1);
subplot(2,1,2);
hold on
% iteration of det fix
iter_stop = [2.0e-3, 2.0e-4];
alpha_iter = 1.0;
Nitermax = 10;
pinfit = pinfit0;
detfix0 = zeros(Npixel*Nslice, Nitermax+1);
for iiter = 1:Nitermax
    % fit splines for each view's projection
    cs = pinsplines(rawdata, fanangles, Npixel, Nslice, Nview);
    % pin matrix
    [A, L, indexrange] = pinleftoptmatrix(cs, Npixel, Nslice, Nview, Npixelpermod, alphaL);
    
    % pin curve fitting
    p0 = reshape(double([cs(:).p]), Nslice, Nview);
    if any(isnan(p0(:)))
        warning('pin out of FOV!');
    end
    
    % fit pin
    [dp, pinfit] = pincurvefit(cs, viewangle, hz, Nslice, Nview, pinfit, coeftofit);
    
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
    fprintf('.');
%     fprintf('%d: %e   %e   %e   %e\n', iiter, norm_detfix, norm_ddf, norminf_detfix, norminf_ddf);
    
    % fix det
    fanangles = fanangles + detfix.*alpha_iter;
    % will be return
    fanfix = fanangles - fanangles0 - detfix.*alpha_iter.*iiter;
    % plot
    if iiter>1
        figure(figid);
        midfix = pinfit(6)*SDD;
        subplot(2,1,1);  imagesc(fanfix'.*SDD - midfix, [-0.15 0.15] - midfix); colormap jet; 
        cb1 = colorbar;
        set(cb1, 'Ticks', [-0.1  0  0.1] - midfix);
        subplot(2,1,2);  plot(mean(detfix, 2));  grid on;
%         ax2 = axis;
%         ax2(1) = 1; ax2(2) = Npixel; 
        ax2 = [1 Npixel -4e-5 4e-5];
        axis(ax2);
        cb2 = colorbar;
        set(cb2, 'Visible', 'off');
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
options = optimoptions('lsqnonlin', 'Display', 'off');
if s(5)     % I know s(5) is the 'rotscale' to scale the viewangles in time based scans 
    [pinfit, ~, ~, exitflag] = lsqnonlin(@(x) pinfitfun2(viewangle, Nslice, hz, x, p0, x0, s), x0, [], [], options);
else
    [pinfit, ~, ~, exitflag] = lsqnonlin(@(x) pinfitfun(viewangle, Nslice, hz, x, p0, x0, s), x0, [], [], options);
end

if exitflag<=0
    if exitflag==0
        warning('please check the pin position!');
    end
    error('pin fitting not converged!');
end
p1 = pinfitfun(viewangle, Nslice, hz, pinfit, 0);
p1 = p1(1:Nslice, 1:Nview);
dp = p1 - p0;
dp = fillmissing(dp, 'constant', 0);

end


function r = pinfitfun(viewangle, Nslice, hz, x, p0, x0, s)

if nargin>5
    x(~s) = x0(~s);
end
x = num2cell(x);

p = pinprojectfun(viewangle, Nslice, hz, x{:});

r = p - p0;
r = fillmissing(r, 'constant', 0);
er = (r - mean(r)).*(Nslice-1);
% dr = r(:,end) - r(:, 1);
r = [r er];

end


function r = pinfitfun2(viewangle, Nslice, hz, x, p0, x0, s)

if nargin>5
    x(~s) = x0(~s);
end
x = num2cell(x);

p = pinprojectfun(viewangle, Nslice, hz, x{:});

r = p - p0;
r = fillmissing(r, 'constant', 0);
er = (r - mean(r)).*(Nslice-1);
dr = r(:,end) - r(:, 1);
r = [r er dr.*size(p, 2)];

end