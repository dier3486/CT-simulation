% run in reconnode_crosstalkcali.m
% try4

% index range
index_range = zeros(2, Nslice, Nview, Nbk/2);
for ii = 1:Nbk/2
    index_range(:,:,:, ii) = reshape(dataflow.(headfields{ii}).index_range, 2, Nslice, Nview);
end

% loop slice
islice = 2;

% Seff
Seff = false(Npixel, Nview*Nbk/2);
for ibk = 1:Nbk/2
    for iview=1:Nview
        viewindex = iview+(ibk-1)*Nview;
        index_ii = index_range(1, islice, iview, ibk) : index_range(2, islice, iview, ibk);
        Seff(index_ii, viewindex) = true;
    end
end

% ipixel = 416;
p = zeros(Npixel, 2);
lamda = 0.1;
options = optimoptions('lsqnonlin','Display','off');
% loop pixel
for ipixel = 2:Npixel-1
    x = zeros(3, Nview*Nbk/2);
    y = zeros(3, Nview*Nbk/2);
    pixelindex = ipixel-1:ipixel+1;
    for ibk = 1:Nbk/2
        ix = ibk*2-1;
        iy = ibk*2;
        viewindex = (1:Nview) + (ibk-1)*Nview;
        x(:, viewindex) = double(dataflow.(datafields{ix})(pixelindex, :));
        y(:, viewindex) = double(dataflow.(datafields{iy})(pixelindex, :));
    end
    s = all(Seff(pixelindex, :), 1);
    if sum(s)<Nview
        continue;
    end
    t0 = [0, 0];
    p(ipixel, :) = lsqnonlin(@(t) crossfit3(t, x, y, s, lamda), t0, [], [], options);
end

