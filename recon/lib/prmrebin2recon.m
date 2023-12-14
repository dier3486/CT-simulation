function recon = prmrebin2recon(recon, rebin)
% to output the recon parameters after rebin prepare
%   recon = prmrebin2recon(recon, rebin);
% call it after rebin prepare.

recon.Nshot = rebin.Nshot;
recon.Npixel = rebin.Nreb;
recon.Nslice = rebin.Nslice;
recon.Nfocal = rebin.Nfocal;
recon.delta_d = rebin.delta_d;
recon.delta_z = rebin.delta_z;
recon.delta_view = rebin.delta_view * rebin.Nfocal;
recon.midchannel = rebin.midU_phi;
recon.SID = rebin.SID;
recon.Nview = rebin.Nview / rebin.Nfocal;
recon.Nviewprot = rebin.Nviewprot / rebin.Nfocal;
recon.gantrytilt = rebin.gantrytilt;

end