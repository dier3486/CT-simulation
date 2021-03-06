function PROTO = protocolcfgsample()
% return a default simulation protocol config, just a sample

% protocal
PROTO.seriesnumber = 2;

% series
PROTO.series{1}.scan = 'Axial';
PROTO.series{1}.collimator = 'All';
PROTO.series{1}.bowtie = 'Body';  % Empty, Body, Head (0 , 1, 2)
PROTO.series{1}.focalspot = 1;
PROTO.series{1}.focalsize = 1;      % 1: small, 2: big
PROTO.series{1}.KV = 120;
PROTO.series{1}.mA = 200;
PROTO.series{1}.mA_air = 30;
PROTO.series{1}.viewperrot = 1440;
PROTO.series{1}.rotationspeed = 1.0;  % sec per rotation
PROTO.series{1}.rotationnumber = 1.0;
PROTO.series{1}.viewnumber = [];
PROTO.series{1}.startangle = 0;
PROTO.series{1}.startcouch = 0;
PROTO.series{1}.shotnumber = 1;
PROTO.series{1}.shotcouchstep = 8.8;
PROTO.series{1}.couchheight = 0;
PROTO.series{1}.couchspeed = 0;
PROTO.series{1}.rawdatastyle = '24bit';

PROTO.series{2}.startcouch = -8.8;
PROTO.series{2}.bowtie = 'Empty';
PROTO.series{2}.shotnumber = 2;

