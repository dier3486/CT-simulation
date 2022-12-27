function interptable = omiga4table(Coeffgamma, Nfourp)
% prepare the table for Sobolev space linearized omiga-4-points interpoloation method
% interptable = omiga4table(Coeffgamma, Nfourp);
% boring code, please look up the omiga4interp.m.

% alpha-beta interp prepare
Nfourp = single(Nfourp);
index_intp = 1:Nfourp+1;
alpha_intp = single(linspace(0, 1, Nfourp+1)');
beta_intp = 1/2-sqrt(1+alpha_intp.*(1-alpha_intp).*4)./2;
fourpoint = [(1+alpha_intp-beta_intp)./2  (alpha_intp+beta_intp)./2  ...
            (Coeffgamma(1)/4)./sqrt(1-alpha_intp.*(1-alpha_intp).*Coeffgamma(2))];

% to return
interptable.Nfourp = Nfourp;
interptable.fourpointindex = index_intp;
interptable.fourpoint = fourpoint;
interptable.convL = single([-1 2 -1]);

end
