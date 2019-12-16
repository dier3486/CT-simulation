function filt = filterdesign(filter, len, d)
% FBP filter
% TBC

% order
order = max(64,2^nextpow2(2*len));

% bandlimited ramp filter 
filt = bandlimitedramp(order);

% frequency axis up to Nyquist
w = 2*pi*(0:size(filt,2)-1)/order; 

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    otherwise
        % handle
        1;
end

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