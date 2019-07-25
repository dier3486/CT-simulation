% test script for configure materials2

addpath(genpath('../'));

material_name = 'water';
material_cfg = loadmaterial(material_name);
samplekeV = (1.0:0.5:90)';
material_def = materialdefine(material_cfg, samplekeV);

figure;
plot(material_def.samplekeV, log([material_def.mu_total, material_def.mu_coh]));