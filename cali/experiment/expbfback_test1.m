% run after bowtieexperiment.m

rsp = load('./testdata/response_1219.mat');
r = rsp.response(:);


% islice = 30;
% iw = 3;

dfit = cell(1, Nw);
dfit(:) = {zeros(Np, Nbow-1)};
dfit_sm = cell(1, Nw);
dfit_sm(:) = {zeros(Npixel, Nslice, Nbow-1)};

mu_1 = SYS.collimation.bowtie{1}.material.mu_total(:);
for iw = 1:Nw
    Dair = log(Pbow{iw}(:,:,1)*(r.*samplekeV(:)));
    Dexp1 = -log(expbow{iw}(:,2)) + log(expbow{iw}(:,1));
    Dexp2 = -log(expbow{iw}(:,3)) + log(expbow{iw}(:,1));
    
    for ipixel = 1:Np
        Pbow1_ip = Pbow{iw}(ipixel, :, 2);
        Pbow2_ip = Pbow{iw}(ipixel, :, 3);
        Dexp1_ip = Dexp1(ipixel)-Dair(ipixel);
        Dexp2_ip = Dexp2(ipixel)-Dair(ipixel);
        
        if isfinite(Dexp1_ip)
            dfit{iw}(ipixel, 1) = fzero(@(x) -log(Pbow1_ip*(exp(-x.*mu_1).*r.*samplekeV(:)))-Dexp1_ip, 0);
        end
        if isfinite(Dexp2_ip)
            dfit{iw}(ipixel, 2) = fzero(@(x) -log(Pbow2_ip*(exp(-x.*mu_1).*r.*samplekeV(:)))-Dexp2_ip, 0);
        end
    end
    dfit{iw} = reshape(dfit{iw}, Npixel, Nslice, Nbow-1);
    for islice = 1:Nslice
        dfit_sm{iw}(:, islice, 1) = smooth(dfit{iw}(:,islice, 1), 0.05, 'rloess');
        dfit_sm{iw}(:, islice, 2) = smooth(dfit{iw}(:,islice, 2), 0.05, 'rloess');
    end
end

% return is dfit, dfit_sm
retdata.bowtiefit = dfit_sm;
retdata.bowtiefit_orig = dfit;
retdata.angle = detangle-pi/2;
% save retdata
