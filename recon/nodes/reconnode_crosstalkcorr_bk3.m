function [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status)
% recon node, crosstalk correction
% [dataflow, prmflow, status] = reconnode_crosstalkcorr(dataflow, prmflow, status);

% parameters to use in prmflow
Nview = prmflow.recon.Nview;
Npixel = prmflow.recon.Npixel;
Nslice = prmflow.recon.Nslice;

% parameters set in pipe
crossprm = prmflow.pipe.(status.nodename);
if isfield(crossprm, 'weight')
    weight = crossprm.weight;
else
    weight = 1.0;
end

% % calibration table
% crscorr = prmflow.corrtable.(status.nodename);
% crsorder = crscorr.order;
% crsval = reshape(crscorr.main, [], crsorder);

% % debug
% Ptest = load('E:\matlab\CT\SINO\PG\calibration\Pcrs_test2.mat');
% crsval = Ptest.pcrs(:);
% crsorder = 1;

% debug2
% load air rate
bhcorr = prmflow.corrtable.Beamharden;
airrate = 2.^(-bhcorr.airrate(:));

Npmod = 16;
alpha = 0.07;
Pmod = ones(Npmod, 1).*alpha;
Pmod(1) = 0;
crsval = repmat(Pmod, Npixel*Nslice/Npmod, 1);
crse1 = ones(Npixel*Nslice, 1).*alpha;
crsorder = 1;

% reshape
dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% to intensity
dataflow.rawdata = 2.^(-dataflow.rawdata);
% correct
% if crsorder == 1
% 	% the correction operator is a tridiagonal matrix [crsval; 1-crsval-crsval_2; crsval_2];
% 	crsval_2 = [crsval(2:end); 0];
% 	% rawfix
%     % diag, A_{i,i} = P_i+P_{i+1}
% 	rawfix = dataflow.rawdata.*(crsval+crsval_2);
%     % diag+1, A_{i,i+1} = -P_{i+1}
% 	rawfix(1:end-1, :) = rawfix(1:end-1, :) - dataflow.rawdata(2:end, :).*crsval(2:end);
%     % diag-1, A_{i,i-1} = -P_i
% 	rawfix(2:end, :) = rawfix(2:end, :) - dataflow.rawdata(1:end-1, :).*crsval(2:end);
% 	% add to rawdata
% 	dataflow.rawdata = dataflow.rawdata + rawfix.*weight;
% else
%     error('Currently thhe crosstalk correction only support 1-order method.');
% end

% for debug
Nps = Npixel*Nslice;
% % style #1: P=D*A*D'
% D = spdiags(repmat([1 -1], Nps, 1), [0 1], Nps, Nps);
% A = spdiags(crsval, 0, Nps, Nps);
% P = D*A*D';
% P(1,1) = -P(1,2);
% style #2: P = (I-Ap)*Ad^(-1)*G, where Ad+Ap=I-D*A*D';
% D = spdiags(repmat([1 -1], Nps, 1), [0 1], Nps, Nps);
% A = spdiags(crsval, 0, Nps, Nps);
% Ap = -D*A*D';
% Ad = spdiags(diag(Ap),0,Nps,Nps);
% Ap = Ap-Ad;
% Ad = Ad+speye(Nps);
% P = (speye(Nps)-Ap)*inv(Ad);
% g = P\ones(Nps,1);
% G = spdiags(g, 0, Nps, Nps);
% P = P*G - speye(Nps);
P = sparsecrossmatrix(crsval);

% Pe1 = sparsecrossmatrix(crse1);
% P = P-Pe1;

% % air rate
alpha = 1.0;
C = spdiags(airrate.^alpha, 0, Nps, Nps);
P = C\(P*C);
rawfix2 = P*double(dataflow.rawdata);
dataflow.rawdata = dataflow.rawdata + rawfix2.*weight;


% min cut
minval = 2^-32;
dataflow.rawdata(dataflow.rawdata<minval) = minval;
% log2
dataflow.rawdata = -log2(dataflow.rawdata);

% status
status.jobdone = true;
status.errorcode = 0;
status.errormsg = [];
end

function C = sparsecrossmatrix(p)

p = p(:);
N = size(p,1);

% diagonal of p
Pdiag = spdiags(p, 0 ,N ,N);
% diff matrix
D = spdiags(repmat([1 -1], N, 1), [0 1], N, N);
% D*P*D'
Ap = -D*Pdiag*D';
Ad = spdiags(diag(Ap), 0, N, N);
Ap = Ap-Ad;
Ad = Ad+speye(N);

C = (speye(N)-Ap)*inv(Ad);
g = C\ones(N,1);
G = spdiags(g, 0, N, N);
C = C*G - speye(N);

end