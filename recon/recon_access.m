function vargout = recon_access(reconxml)
% recon & cali governing function
% output = recon_access(reconxml);

if ischar(reconxml)
    % try to read recon xml file
    reconxml = readcfgfile(reconxml);
end

if ~iscell(reconxml.recon)
    reconxml.recon = {reconxml.recon};
end



end
