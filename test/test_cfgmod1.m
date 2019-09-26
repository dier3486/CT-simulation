% test to define MOD
clear

phantom_cfg = phantomcfgsample();

root.test_phantom = phantom_cfg;
xmlcfgfile = './test phantom.xml';
struct2xml(root, xmlcfgfile);

% 
% jsontxt = jsonwrite(phantom_cfg);

p1 = phantomconfigure(xmlcfgfile);

samplekeV = [50, 100];
p2 = materialconfigure(p1, samplekeV);


