% the color palette of Z-map

status.nodename = 'Materialzmap';

% reset prm
prmflow.pipe.Materialzmap = struct();
prmflow.pipe.Materialzmap.ImageWindow = [50 5050];
prmflow.pipe.Materialzmap.LightLorenz = 0;
prmflow.pipe.Materialzmap.LightRefWindow = 2000;
prmflow.pipe.Materialzmap.LightMin = 0.1;
prmflow.pipe.Materialzmap.Saturation = 0.75;


prmflow.pipe.Materialzmap.showpalette = true;

% replay prepare
[dataflow, prmflow, status] = reconnode_materialzmapprepare(dataflow, prmflow, status);

% Z-map
[dataflow, prmflow, status] = reconnode_materialzmap(dataflow, prmflow, status);

% show images
cttool(dataflow.imageColor);