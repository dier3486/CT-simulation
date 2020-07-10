function img = imagereorder(img, Nrow, reorderflag)
% reorder and ' the images after BP
% img = imagereorder(img, Nrow, reorderflag);
% normally, img is dataflow.image{iseries};
% Nrow is prmflow.recon.Nslice;
% reorderflag is protocol.couchdirection<0 & Axial

if reorderflag
    % do nothing
    return;
end

Nimg = size(img, 3);
index = flipud(reshape(1:Nimg, Nrow, []));
img = img(:,:,index);

end