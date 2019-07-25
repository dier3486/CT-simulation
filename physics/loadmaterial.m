function material_cfg = loadmaterial(material_file)
% load material configure file
% material_cfg = loadmaterial(material_file) or loadmaterial(material_name)

if exist(material_file, 'file')
    material_cfg = jsonread(material_file);
elseif exist([material_file, '.cfg'], 'file')
    material_cfg = jsonread([material_file, '.cfg']);
else
    error(['Can not find .cfg file for material: ', material_file]);
end
