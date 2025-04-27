function A = putfieldsinGPU(A, fields)
% put fields of stuct A in GPU 

for ii = 1:length(fields)
    if isfield(A, fields{ii})
        A.(fields{ii}) = everything2single(A.(fields{ii}), 'any', 'gpuArray');
    end

end