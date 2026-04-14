function [datainfo, CUdatainfo] = pooldatainfo(poolsize, pooldata, fieldkey, Index)
% to return the pooldatainfo used in CUpipeline.h, it is a parepare of
% data copy.

datainfo = struct();
datainfo.poolsize = poolsize;
if isfield(pooldata, [fieldkey '_databasis'])
    datainfo.databasis = pooldata.([fieldkey '_databasis']);
else
    datainfo.databasis = 1;
end
if isfield(pooldata, [fieldkey '_datasize'])
    datainfo.datasize = pooldata.([fieldkey '_datasize']);
else
    datainfo.datasize = 1;
end
if isfield(pooldata, [fieldkey '_dataclass'])
    datainfo.dataclass = pooldata.([fieldkey '_dataclass']);
else
    datainfo.dataclass = 'single';
end
datainfo.classsize = classsize(datainfo.dataclass);
if nargin > 3
    datainfo.index = Index;
else
    datainfo.index = [1 0];
end

if nargout > 1 && libisloaded('CUpipeline')
    datainfo2 = datainfo;
    datainfo2.datasize = datainfo.classsize * datainfo.datasize * datainfo.databasis;
    datainfo2.index = datainfo.index - 1;
    CUdatainfo = libstruct('pooldatainfo', datainfo2);
else
    CUdatainfo = [];
end

end