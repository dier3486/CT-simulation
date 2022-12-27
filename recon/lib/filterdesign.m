function filt = filterdesign(filter, len, d, a, varargin)
% FBP filter design
% filt = filterdesign(filter, len, d, a)

if nargin<4
    a = 1.0;
end

% order
order = max(64,2^nextpow2(len*1.2));

% bandlimited ramp filter 
filt = bandlimitedramp(order);

% frequency axis up to Nyquist
w = 2*pi*(0:size(filt,2)-1)/order; 

if isempty(filter)
    % call ram-lak as default filter
    filter = 'ram-lak';
end

switch lower(filter)
    case 'none'
        % no filter
        filt = ones(order, 1);
        return
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d*a))./(w(2:end)/(2*d*a)));
        filt(filt<0) = 0;
%         filt(2:end) = filt(2:end) .* (sin(w(2:end)/(1.2*d))./(w(2:end)/(1.2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d*a));
        filt(filt<0) = 0;
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/(d*a)));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./(d*a))) / 2;
    case 'fermi'
        alpha = log(19)/varargin{1};
        filt(12:end) = filt(12:end) ./(1+exp((w(12:end)-(d*a)*pi/2).*alpha));
        1;
    case 'myfilter'
        % TBC
    otherwise
        % handle
        % TBC
        1;
end

filt = filt./d;
filt(w>pi*(d*a)) = 0; 
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter

end

function filt = bandlimitedramp(order)

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

end