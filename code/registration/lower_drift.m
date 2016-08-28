function Iout = lower_drift(I,sigma,debug)
    sz = ceil((sigma*6 + 1)/2)*2 - 1;
    f = fspecial('gaussian',[sz sz], sigma);
    tmp = imfilter(I,f);
    Iout = mat2gray(I - tmp);
    if debug
        figure,imshowpair(Iout,f,'montage');
        figure
        subaxis(2,2,1); imagesc(I); axis image; colormap gray
        subaxis(2,2,2); imagesc(tmp); axis image; colormap gray
        subaxis(2,2,3); imagesc(Iout); axis image; colormap gray
        %subaxis(2,2,4); imagesc(adapthisteq(tmp1)); axis image; colormap gray
        subaxis(2,2,4); imagesc(edge((Iout),'canny',[])); axis image; colormap gray
    end
end