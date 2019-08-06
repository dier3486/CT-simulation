function MOD_cfg = loadmodconfig(cfgfile)
% load MOD configure file, MOD is the objects to scan in CT 
% MOD_cfg = modconfig(cfgfile)



return

function object = default_object()
object.name = 'default';
object.parent = 0;
object.type = '';
object.O = [0, 0, 0];
object.vector = eye(3);
object.invV = eye(3);
object.shape_adv = struct();
object.material = '';



return
