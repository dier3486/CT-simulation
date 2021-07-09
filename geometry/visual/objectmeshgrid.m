function [X, Y, Z, C] = objectmeshgrid(object, varargin)
% get the XYZ menshgrid of an object to plot(mesh)
% [X, Y, Z] = objectmeshgrid(object)
% or [X, Y, Z, C] = objectmeshgrid(object, nshell, nangle, ncover);

switch lower(object.type)
    case {'sphere', 'ellipsoid'}
        % sphere
        [nshell, nangle] = cleaninputarg({20, 24}, varargin{:});
        thetagrid = linspace(0, pi*2, nangle+1);
        phigrid = linspace(0, pi, nshell);
        X = sin(phigrid(:))*cos(thetagrid);
        Y = sin(phigrid(:))*sin(thetagrid);
        Z = repmat(cos(phigrid(:)), 1, nangle+1);   
    case {'cylinder'}
        % cylinder
        [nshell, nangle, ncover] = cleaninputarg({10, 24, 5}, varargin{:});
        thetagrid = linspace(0, pi*2, nangle+1);
        zshell = linspace(-1, 1, nshell);
        zgrid = [-ones(1, ncover) zshell ones(1, ncover)];
        rhocov = linspace(0, 1, ncover+1);
        rhogrid = [rhocov(1:end-1) ones(1, nshell) rhocov(end-1:-1:1)];
%         [theta, rho] = meshgrid(thetagrid, rhogrid);
        X = rhogrid(:)*cos(thetagrid);
        Y = rhogrid(:)*sin(thetagrid);
        Z = repmat(zgrid(:), 1, nangle+1);
    case 'tube'
        % cylinder without cover
        [nshell, nangle] = cleaninputarg({12, 24}, varargin{:});
        thetagrid = linspace(0, pi*2, nangle+1);
        zgrid = linspace(-1, 1, nshell);
        X = repmat(cos(thetagrid), nshell, 1);
        Y = repmat(sin(thetagrid), nshell, 1);
        Z = repmat(zgrid(:), 1, nangle+1);
    case 'cone'
        % cone
        if isfield(object, 'anglerange')
            anglerange = object.anglerange;
        else
            anglerange = [0 pi*2];
        end
        nangle = ceil(abs(anglerange(2)-anglerange(1))/pi*12);
        [nshell, nangle] = cleaninputarg({12, nangle}, varargin{:});
        thetagrid = linspace(anglerange(1), anglerange(2), nangle+1);
        zgrid = linspace(0, 1, nshell);
        X = zgrid(:)*cos(thetagrid);
        Y = zgrid(:)*sin(thetagrid);
        Z = repmat(zgrid(:), 1, nangle+1);
    otherwise
        warning('Unknown object type %s!', object.type);
        X = []; Y = []; Z = [];
end

C = Z;
gridsize = size(Z);
[X, Y, Z] = tac(reshape([X(:) Y(:) Z(:)]*object.vector+object.O, [gridsize 3]), [1 2]);

end

