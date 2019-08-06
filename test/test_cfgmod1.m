% test to define MOD
clear

MOD_cfg.name = 'test1';
MOD_cfg.rootnode = 0;
MOD_cfg.nodetree = [0, 1, 1];

% default_obj.name = 'default';
% defualt_obj.parent = 0;
default_obj.type = ' ';
default_obj.O = [0, 0, 0];
default_obj.vector = eye(3);
default_obj.invV = [];
defualt_obj.volume = [];
default_obj.shape_adv = struct();
default_obj.material = 'xxx';

MOD_cfg.object = cell(3,1);
MOD_cfg.object(:) = {default_obj};
% MOD_cfg.object(1:3) = default_obj;

m1.MOD_cfg = MOD_cfg;
% xml1=struct2xml(m1)