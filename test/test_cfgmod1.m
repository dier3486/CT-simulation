% test to define MOD
clear

phantom_cfg = phantomcfgsample();

xml.phantom_cfg = phantom_cfg;
xmltxt=struct2xml(xml);

jsontxt = jsonwrite(phantom_cfg);


