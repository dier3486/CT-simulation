function BC = TVboundarycond(TVopts)
BC = zeros(3, 2);
% BC(dim, left/right)

bcfields = {'BoudaryCondX', 'BoudaryCondY', 'BoudaryCondZ'};

for ibc = 1:3
    switch(TVopts.(bcfields{ibc}))
        case "fix"
            BC(ibc, :) = 0;
        case "clamp"
            BC(ibc, :) = 1;
        case "linextrap"
            BC(ibc, :) = 2;
        case "mirror"
            BC(ibc, :) = 3;
        case "oddextrap"
            BC(ibc, :) = 4;
        otherwise
            error(0);
    end
end

% The left/right boundary conditions are same now, but
if TVopts.LeftBoundary
    BC(TVopts.LeftDim, 1) = 0;
end


end