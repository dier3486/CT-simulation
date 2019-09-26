% test script for bin file io
addpath(genpath('../'));

S1 = struct();

S1(1).a = 1.0;
S1(2).a = pi;

[S1(:).b] = deal(struct());
S1(1).b(1).v = single(1.0);
S1(1).b(1).c = 'abc';
S1(1).b(2).v = single(2.0);
S1(1).b(2).c = 'cba';

S1(2).b(1).v = single(6.0);
S1(2).b(1).c = 'ABC';
S1(2).b(2).v = single(8.0);
S1(2).b(2).c = 'CBA';

S1(1).c = 'ab';
S1(2).c = 'cd';

S1(1).d = uint8([1 2 3 4]);
S1(2).d = uint8([3 5 7 9]);

S1(1).l = logical([1 0 1 0]);
S1(2).l = logical([1 1 1 1]);

cfg1 = structbincfg(S1);

[data1, cfg1_p] = packstruct(S1, cfg1);

fid = fopen('./testdata/data1.bin', 'w');
fwrite(fid, data1, 'uint8');
fclose(fid);

root = struct();
root.S = cfg1;
xmlfile = './testdata/data1.xml';
struct2xml(root, './testdata/data1.xml');
jsonfile = './testdata/data1.json';
jsonwrite(cfg1, jsonfile);

S2 = sparsepack(data1, cfg1);