hfg = figure;
himg = image();
A = himg.CData;
himg.delete;
close(hfg);

% A = zeros(64, 64);
% A(50, 32) = 1;

Na = size(A, 1);
Xa = (-(Na-1)/2:(Na-1)/2);
Ya = Xa;

Nv = 180;
Vb = (0:Nv-1).*pi/Nv - pi/2;

Nb = ceil(Na*sqrt(2)/2)*2+1;
h = 1;
Rb = (-(Nb-1)/2 : (Nb-1)/2).*h;
[ndRb, ndVb] = ndgrid(Rb, Vb);

Xb = ndRb.*cos(ndVb);
Yb = ndRb.*sin(ndVb);

B = interp2(Xa, Ya, A', Xb, Yb);

Bd = B(2:end-1, :) - B(1:end-2, :)./2 - B(3:end, :)./2;

[ndXa, ndYa] = ndgrid(Xa, Ya);

Va = atan(ndYa./ndXa);
Ra = sqrt(ndYa.^2 + ndXa.^2).*sign(ndXa);

% A1 = interp2([Vb pi/2], Rb, [B B(:,1)], Va, Ra);

B1 = mean(Bd, 2, 'omitnan');
B1 = (B1 + flipud(B1))./2;
B1 = repmat(B1, 1, Nv+1);
C1 = interp2([Vb pi/2], Rb(2:end-1), B1, Va, Ra);