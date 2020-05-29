function img = loaddcmimage(imgpath, dcmext)

if nargin<2
    dcmext = '.dcm';
end

files = dir(fullfile(imgpath, [ '*' dcmext]));
[~, sortidx] = natsortfiles({files.name});
files = files(sortidx);
Nf = size(files(:), 1);

img = [];
for ii = 1:Nf
    img_ii = dicomread(fullfile(files(ii).folder, files(ii).name));
    if ii == 1
        img = zeros([size(img_ii) Nf]);
    end
    img(:,:,ii) = img_ii;
end

end