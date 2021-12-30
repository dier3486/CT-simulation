function imgout = enhancemix(img1, img2, wlevel, mixwidth, alpha, beta)

t = img1-wlevel;
t(t>mixwidth/2) = mixwidth/2;
t(t<-mixwidth/2) = -mixwidth/2;

a = mixwidth^2*3/4;
b = beta*2/mixwidth^3;
fmix = t.*(a-t.^2).*b + alpha;

imgout = img1.*fmix + img2.*(1-fmix);

end