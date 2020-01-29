function Dmu = projectinphantom(focalposition, detectorposition, phantom, samplekeV, viewangle, couch)
% project in phantom
% Dmu = projectinphantom(focalposition, detectorposition, phantom, samplekeV, viewangle, couch);

% default viewangle(s) and couch
if nargin<5
    viewangle = 0;
    couch = 0;
end
if nargin<6
    couch = zeros(size(viewangle));
end

% size
Np = size(detectorposition, 1);
Nview = length(viewangle(:));
Nsample = length(samplekeV(:));
Nfocal = size(focalposition, 1);
viewangle = reshape(viewangle, Nfocal, []);
% couch = reshape(couch, Nfocal, [], 3);

% ini Dmu
Dmu = zeros(Np*Nview, Nsample);

if ~isempty(phantom)
    for iobj = 1:phantom.Nobject
        % loop the objects
        parentobj = phantom.object_tree(iobj);
        object_i = phantom.object{iobj};
        % mu
        mu_i = interp1(object_i.material.samplekeV, object_i.material.mu_total, samplekeV);
        % mu fix
        if parentobj>0
            mu_parent = interp1(phantom.object{parentobj}.material.samplekeV, ...
                phantom.object{parentobj}.material.mu_total, samplekeV);
            mu_i = mu_i - mu_parent;
        end
        % ini D_i
        D_i = zeros(Np, Nview);
        % L = zeros(Np, Nfocal);
        % I know the L has been done in flewoverbowtie.m
        % fly focal
        for ifocal = 1:Nfocal
            % geometry projection in object
            [D_i(:, ifocal:Nfocal:end), ~] = intersection(focalposition(ifocal, :), ...
                detectorposition, object_i, 'views-ray', viewangle(ifocal, :), couch(ifocal:Nfocal:end, :)); 
        end
        Dmu = Dmu + D_i(:)*mu_i;
    end
else
    % donothing
    1;
end

end
