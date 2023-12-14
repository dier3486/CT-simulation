function [FPchannelpos, effNp, indexstart] = fpprepare(recon, viewangle)
% a subfunction of forward projection prepare 
%   [FPchannelpos, effNp, indexstart] = fpprepare(recon, viewangle);
% or, [FPchannelpos, effNp] = fpprepare(recon);

% if isfield(recon, 'upsampled')
%     upsampled = recon.upsampled;
% else
%     upsampled = false;
% end

% FPchannelpos
FPchannelpos = ((1:recon.Npixel)'-recon.midchannel).*recon.delta_d;
effNp = ceil(recon.effFOV/recon.delta_d) + 2;

% viewangle shall be dataflow.head.viewangle, do not use recon.viewangle! 
if nargin>1
    eta_C = recon.center(1).*sin(viewangle) - recon.center(2).*cos(viewangle);
    indexstart_p = floor(recon.midchannel + (-recon.effFOV/2 + eta_C)./recon.delta_d);
    indexstart_n = floor(recon.midchannel*2 - indexstart_p);
    indexstart = [indexstart_p(:)  indexstart_n(:)];
end


end