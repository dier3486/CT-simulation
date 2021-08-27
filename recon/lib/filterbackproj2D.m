function [img,H] = filterbackproj2D(p, projprm, filter)
% copy from IRADON Inverse Radon transform.
% 

%   Copyright 1993-2013 The MathWorks, Inc.

%   References:
%      A. C. Kak, Malcolm Slaney, "Principles of Computerized Tomographic
%      Imaging", IEEE Press 1988.

% modified by Dier Z.
% It is a test code, do not call it in recon nodes.

theta = projprm.viewangle;
ctrIdx = projprm.midchannel;
dscale = projprm.delta_d;
hond = projprm.h/projprm.delta_d;
if isfield(projprm, 'imagesize')
    N = projprm.imagesize;
else
    N = 512;
end
if isfield(projprm, 'interp')
    interp = projprm.interp;
else
    interp = 'linear';
end
if isfield(projprm, 'Nslice')
    Nslice = projprm.Nslice;
else
    Nslice = 1;
end
if isfield(projprm, 'Nview')
    Nview = projprm.Nview;
else
    Nview = size(p, 2);
end
if isfield(projprm, 'maxR')
    maxR = projprm.maxR;
else
    maxR = inf;
end
if nargin<3
    filter = 'ram-lak';    % The ramp filter is the default
end

p = reshape(p, projprm.Npixel, []);
[p,H] = filterProjections(p, filter, dscale);


% len = size(p,1);
% ctrIdx = ceil(len/2);     % index of the center of the projections

% Zero pad the projections to size 1+2*ceil(N/sqrt(2)) if this
% quantity is greater than the length of the projections
imgDiag = 2*ceil(N/sqrt(2))+1;  % largest distance through image.
if size(p,1) < imgDiag
    rz = imgDiag - size(p,1);  % how many rows of zeros
    p = [zeros(ceil(rz/2),size(p,2)); p; zeros(floor(rz/2),size(p,2))];
    ctrIdx = ctrIdx+ceil(rz/2);
end

p = reshape(p, [], Nslice, Nview);
if isa(p, 'gpuArray')
    img = backproj2D_GPU(p, theta, ctrIdx, hond, N, [0 0], maxR);
else
    img = backproj2D_2(p, theta, ctrIdx, hond, N, interp);
end

if isfield(projprm, 'fillmiss')
    img = fillmissing(img, 'const', projprm.fillmiss);
end

return

%======================================================================
function [p,H] = filterProjections(p_in, filter, d)

p = p_in;

% Design the filter
len = size(p,1);
if ischar(filter)
    H = designFilter(filter, len, d);
else
    H = filter;
end

if strcmpi(filter, 'none')
    return;
end

p(length(H),1)=0;  % Zero pad projections

% In the code below, I continuously reuse the array p so as to
% save memory.  This makes it harder to read, but the comments
% explain what is going on.

p = fft(p);               % p holds fft of projections

p = bsxfun(@times, p, H); % faster than for-loop

p = ifft(p,'symmetric');  % p is the filtered projections
% p = ifft(p);

p(len+1:end,:) = [];      % Truncate the filtered projections
%----------------------------------------------------------------------

%======================================================================
function filt = designFilter(filter, len, d)
% Returns the Fourier Transform of the filter which will be
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


order = max(64,2^nextpow2(2*len));

if strcmpi(filter, 'none')
    filt = ones(1, order);
    return;
end

% First create a bandlimited ramp filter (Eqn. 61 Chapter 3, Kak and
% Slaney) - go up to the next highest power of 2.

n = 0:(order/2); % 'order' is always even. 
filtImpResp = zeros(1,(order/2)+1); % 'filtImpResp' is the bandlimited ramp's impulse response (values for even n are 0)
filtImpResp(1) = 1/4; % Set the DC term 
filtImpResp(2:2:end) = -1./((pi*n(2:2:end)).^2); % Set the values for odd n
filtImpResp = [filtImpResp filtImpResp(end-1:-1:2)]; 
filt = 2*real(fft(filtImpResp)); 
filt = filt(1:(order/2)+1);

w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist

switch filter
    case 'ram-lak'
        % Do nothing
    case 'shepp-logan'
        % be careful not to divide by 0:
        filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
    case 'cosine'
        filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
    case 'hamming'
        filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
    case 'hann'
        filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
    otherwise
        error(message('images:iradon:invalidFilter'))
end
filt = filt./d;
filt(w>pi*d) = 0;                      % Crop the frequency response
filt(filt<0) = 0;
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
%----------------------------------------------------------------------


%======================================================================
function [p,theta,ctrIdx,filter,d,interp,N] = parse_inputs(varargin)
%  Parse the input arguments and retun things
%
%  Inputs:   varargin -   Cell array containing all of the actual inputs
%
%  Outputs:  p        -   Projection data
%            theta    -   the angles at which the projections were taken
%            ctrIdx   -   % index of the center of the projections
%            filter   -   string specifying filter or the actual filter
%            d        -   a scalar specifying normalized freq. at which to crop
%                         the frequency response of the filter
%            interp   -   the type of interpolation to use
%            N        -   The size of the reconstructed image

p     = varargin{1};
theta = varargin{2};

validateattributes(p    ,{'numeric','logical'},{'real','2d'},mfilename,'R'    ,1);
validateattributes(theta,{'numeric','logical'},{'real'}     ,mfilename,'theta',2);

% % Convert to radians % delete
% theta = pi*theta/180;

% Default values
N = 0;                 % Size of the reconstructed image
d = 1;                 % Defaults to no cropping of filters frequency response
filter = 'ram-lak';    % The ramp filter is the default
interp = 'linear';     % default interpolation is linear

interp_strings = {'nearest neighbor', 'linear', 'spline', 'pchip', 'cubic', 'v5cubic'};
filter_strings = {'ram-lak','shepp-logan','cosine','hamming', 'hann', 'none'};
string_args    = [interp_strings filter_strings];

for i=3:nargin
    arg = varargin{i};
    if ischar(arg)
        str = validatestring(arg,string_args,mfilename,'interpolation or filter');
        idx = find(strcmp(str,string_args),1,'first');
        if idx <= numel(interp_strings)
            % interpolation method
            interp = string_args{idx};
        else %if (idx > numel(interp_strings)) && (idx <= numel(string_args))
            % filter type
            filter = string_args{idx};
        end
    elseif numel(arg)==1
        if arg <=1
            % frequency scale
            validateattributes(arg,{'numeric','logical'},...
                {'positive','real','nonsparse'},mfilename,'frequency_scaling');
            d = arg;
        else
            % output size
            validateattributes(arg,{'numeric'},{'real','finite','integer'},...
                mfilename,'output_size');
            N = arg;
        end
    else
        error(message('images:iradon:invalidInputParameters'))
    end
end

% If 'cubic' interpolation is specified, convert it to 'pchip'. 
if strcmp(interp,'cubic')
    interp = 'pchip';
end 

% If the user didn't specify the size of the reconstruction, so
% deduce it from the length of projections
if N==0
    N = 2*floor( size(p,1)/(2*sqrt(2)) );  % This doesn't always jive with RADON
end

% for empty theta, choose an intelligent default delta-theta
if isempty(theta)
    theta = pi / size(p,2);
end

% If the user passed in delta-theta, build the vector of theta values
if numel(theta)==1
    theta = (0:(size(p,2)-1))* theta;
end

if length(theta) ~= size(p,2)
    error(message('images:iradon:thetaNotMatchingProjectionNumber'))
end
%----------------------------------------------------------------------
