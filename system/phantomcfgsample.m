function phantom_cfg = phantomcfgsample()
% return a sample CT phantom

phantom_cfg.Nobject = 3;
phantom_cfg.object_root = 0;
phantom_cfg.object_tree = [0 1 1];

phantom_cfg.object{1} = struct();
phantom_cfg.object{1}.type = 'Cylinder';
phantom_cfg.object{1}.O = [0, 0, 0];
v1 = [160  0  0; 0 160 0; 0 0 30];
phantom_cfg.object{1}.vector = v1;
% phantom_cfg.object{1}.volume = defaultvolume(v1, phantom_cfg.object{1}.type);
% phantom_cfg.object{1}.material = loadmaterial('water');
phantom_cfg.object{1}.material = 'Water';

phantom_cfg.object{2} = struct();
phantom_cfg.object{2}.type = 'Sphere';
phantom_cfg.object{2}.O = [120, 0, 0];
v2 = [20  0  0; 0 20 0; 0 0 20];
phantom_cfg.object{2}.vector = v2;
% phantom_cfg.object{2}.volume = defaultvolume(v2, phantom_cfg.object{2}.type);
% phantom_cfg.object{2}.material = loadmaterial('vacuum');
phantom_cfg.object{2}.material = 'Vacuum';

phantom_cfg.object{3} = struct();
phantom_cfg.object{3}.type = 'Sphere';
phantom_cfg.object{3}.O = [0, 80, 6];
v3 = [40  0  0; 0 40 0; 0 0 20];
phantom_cfg.object{3}.vector = v3;
% phantom_cfg.object{3}.volume = defaultvolume(v3, phantom_cfg.object{2}.type);
% phantom_cfg.object{3}.material = loadmaterial('water1100');
phantom_cfg.object{3}.material = 'Water1100';
