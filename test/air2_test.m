% test 2.0 ari correction

% % parameters to use in prmflow
% Npixel = prmflow.recon.Npixel;
% Nslice = prmflow.recon.Nslice;
% Nview = prmflow.recon.Nview;
% Nfocal = prmflow.system.Nfocal;
% 
% % calibration table
% aircorr = prmflow.corrtable.Air;
% % parameters in corr
% Nsect = single(aircorr.Nsection);
% refpixel = single(aircorr.refpixel);
% Nref = aircorr.refnumber;
% 
% % angles of the air table
% sectangle = (pi*2/Nsect);
% % airangle = (-1:Nsect).*(pi*2/Nsect);
% % airmain & airref
% aircorr.main = reshape(aircorr.main, [], Nsect);
% airmain = [aircorr.main aircorr.main(:,1)];
% aircorr.reference = reshape(aircorr.reference, [], Nsect);
% airref = [aircorr.reference aircorr.reference(:,1)];
% 
% % interp index and weight
% retangle = mod(dataflow.rawhead.viewangle - aircorr.firstangle, pi*2);
% intp_index = floor(retangle./sectangle);
% intp_alpha = retangle./sectangle - intp_index;
% intp_index = intp_index + 1;
% 
% % reshape
% dataflow.rawdata = reshape(dataflow.rawdata, Npixel*Nslice, Nview);
% 
% % corr rawdata with air
% for ifocal = 1:Nfocal
%     viewindex = ifocal:Nfocal:Nview;
%     airindex = (1:Npixel*Nslice) + Npixel*Nslice*(ifocal-1);
%     dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index).*(1-intp_alpha(viewindex));
%     dataflow.rawdata(:, viewindex) = dataflow.rawdata(:, viewindex) - airmain(airindex, intp_index+1).*intp_alpha(viewindex);
% end
% 
% % skip the edge slices
% if Nslice>2
%     index_slice = 2:Nslice-1;
%     Nrefsl = Nslice-2;
% else
%     index_slice = 1:Nslice;
%     Nrefsl = Nslice;
% end
% 
% dataflow.rawdata = reshape(dataflow.rawdata, Npixel, Nslice, []);
% 
% % referenece data
% ref1 = reshape(dataflow.rawdata(1:refpixel, index_slice, :), [], Nview);
% ref2 = reshape(dataflow.rawdata(Npixel-refpixel+1:Npixel, index_slice, :), [], Nview);
% 
% % reference error
% ref1_err = reshape(ref1-mean(ref1, 1), Nrefsl, refpixel, Nview);
% ref2_err = reshape(ref2-mean(ref2, 1), Nrefsl, refpixel, Nview);
% % SVD 
% referr = zeros(2, Nview);
% for iview = 1:Nview
%     s1 = svd(ref1_err(:, :, iview));
%     referr(1, iview) = s1(1);
%     s2 = svd(ref2_err(:, :, iview));
%     referr(2, iview) = s2(1);
% end
% % norm
% referr = referr./sqrt(Nrefsl*refpixel);
% 
% % rawref
% ref1 = mean(ref1, 1);
% ref2 = mean(ref2, 1);
% % rawref = [ref1; ref2];

% ref block
block_cut = 0.05;
m_blk = 4;
blk1 = conv(ref1>block_cut, ones(1, 2*m_blk+1));
blk2 = conv(ref2>block_cut, ones(1, 2*m_blk+1));
blk1 = blk1(m+1:end-m)>0;
blk2 = blk2(m+1:end-m)>0;

idx_both = blk1 & blk2;

ref = zeros(1, Nview);
ref(idx_both) = (ref1(idx_both) + ref2(idx_both))./2;

idx_1 = blk1 & ~blk2;
idx_2 = ~blk1 & blk2;

