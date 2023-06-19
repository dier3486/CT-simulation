function sampledata = spectrumresample(spectsample, spectdata, samplekeV)
% resampling the spectrum data by the samplekeV
%   sampledata = spectrumresample(spectsample, spectdata, samplekeV);
% typically
%   source.spectrum{ii} = spectrumresample(spectdata(:,1), spectdata(:,2), samplekeV);
% It is not interplation.

dx = mean(diff(spectsample));
dy = mean(diff(samplekeV));

Nx = length(spectsample);
Ny = length(samplekeV);

spectdata = [0; spectdata(:); 0];
t0 = (samplekeV(1) - spectsample(1) - dy + dx)/dx;
d = dy/dx;
sampledata = samplekeV.*0;

for ii = 1:Ny
    tm = t0 + (ii-1/2)*d + 1/2;
    tn = t0 + (ii+1/2)*d - 1/2;

    m = max(floor(tm), 0);
    n = min(ceil(tn), Nx+1);
    if n<0
        % out of lower boundary
        continue;
    elseif m>Nx+1
        % out of upper boundary
        break;
    end
    sampledata(ii) = sum(spectdata(m+1:n+1)) - spectdata(m+1)*(tm-m) + spectdata(n+1)*(tn-n);
end

end