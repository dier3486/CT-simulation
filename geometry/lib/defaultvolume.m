function vol = defaultvolume(V, type)
% retrun volume of an object
% vol = defaultvolume(object.vector, object.type);

switch lower(type)
    case {'cube', 'paralleldpiped'}
        vol = det(V);
    case {'sphere', 'ellipsoid'}
        vol = det(V)*(pi*3/4);
    case {'cylinder'}
        vol = det(V)*(pi*2);
    otherwise
        vol = 0;
end

return