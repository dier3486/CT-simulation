function material_def = materialdefine(material_cfg, samplekeV, elementspath)
% to define a material 
% material_def = materialdefine(material_cfg, samplekeV, elementspath)
% where the material_cfg can read by loadmaterial(material_file)

if nargin < 2
    samplekeV = [];
end
if nargin < 3
    elementspath = '';
end

% copy
material_def.name = material_cfg.name;
material_def.density = material_cfg.density;
material_def.elements = fieldnames(material_cfg.elements);
material_def.Nelem = size(material_def.elements,1);
material_def.elemweight = zeros(material_def.Nelem, 1);
material_def.elemdata = {};
% load elements data
for i_elem = 1:material_def.Nelem
    material_def.elemprm(i_elem) = jsonread([elementspath, material_def.elements{i_elem}, '.txt']);
    material_def.elemdata{i_elem} = csvread([elementspath, material_def.elements{i_elem}, '.csv'], 1, 0);
    switch material_cfg.method
        case 'mol'
            material_def.elemweight(i_elem) = material_cfg.elements.(material_def.elements{i_elem}) * ...
                material_def.elemprm(i_elem).weight;
        case 'weight'
            material_def.elemweight(i_elem) = material_cfg.elements.(material_def.elements{i_elem});
        otherwise
            % error
            error(['Unknown material defination method: ''', material_cfg.method, '''']);
    end
end
material_def.elemweight = material_def.elemweight./sum(material_def.elemweight);

% interpolate the data to samplekeV
samplekeV = samplekeV(:);
material_def.samplekeV = samplekeV;
material_def.mu_total = zeros(size(samplekeV));
material_def.mu_coh = zeros(size(samplekeV));
% I know
index_keV = 1;
index_total = 6;
index_coh = 5;
% loop the elements
for i_elem = 1:material_def.Nelem
    edges = [-inf, material_def.elemprm(i_elem).E_edges(:)', inf];
    % loop the E edges
    for i_edge = 1:length(edges)-1
        % interpolation
        s_sample = samplekeV>=edges(i_edge) & samplekeV<edges(i_edge+1);
        if sum(s_sample)==0
            continue
        end
        s_data = material_def.elemdata{i_elem}(:, index_keV) >= edges(i_edge) ...
                 & material_def.elemdata{i_elem}(:, index_keV) < edges(i_edge+1);
        mutotal_ielem = interp1(material_def.elemdata{i_elem}(s_data, index_keV), ...
            material_def.elemdata{i_elem}(s_data, index_total), samplekeV(s_sample), 'spline');
        mucoh_ielem = interp1(material_def.elemdata{i_elem}(s_data, index_keV), ...
            material_def.elemdata{i_elem}(s_data, index_coh), samplekeV(s_sample), 'spline');
        material_def.mu_total(s_sample) = material_def.mu_total(s_sample) + ...
            mutotal_ielem.*material_def.elemweight(i_elem);
        material_def.mu_coh(s_sample) = material_def.mu_coh(s_sample) + ...
            mucoh_ielem.*material_def.elemweight(i_elem);
    end
end

