function object = fillupobject(object)
% fillup an object

if ~isstruct(object)
    return
end

if isfield(object, 'invV')
    object.invV = inv(object.vector);
elseif isempty(object.invV)
    object.invV = inv(object.vector);
end

if isfield(object, 'volume')
    object.volume = defaultvolume(object.vector, object.type);
elseif isempty(object.volume)
    object.volume = defaultvolume(object.vector, object.type);
end

return

