function Imageout = GetBoneEdge_test1(ImgIn,minValue)

  Imageout =  ImgIn;


   w=fspecial('log',15,0.45);
   w = -w;
   img =ImgIn;
  img(img < minValue) = 0;
    %  imtool(img);
   img(img >= minValue) = 1;
     %  imtool(img);
   img = filter2(w,img);
      % imtool(img);
   img(img < 0) = 0;
    img = img+1;
%    imtool(img);
   Imageout = (Imageout-minValue).*img + minValue;
         
end
