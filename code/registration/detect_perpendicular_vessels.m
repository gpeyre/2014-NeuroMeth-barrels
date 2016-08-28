function [centroids,corrs,radius] = detect_perpendicular_vessels(I1,sigma,thresh,blackAndWhite)
    [N,M] = size(I1);
    nsigma = size(sigma,2);
    cdetect = zeros(N,M,nsigma);
    %figure
    for i = 1:nsigma
        % Detect white bumps
        sz = ceil((sigma(i)*6 + 1)/2)*2 - 1;% size of the gaussian template
        f = fspecial('gaussian',[sz sz], sigma(i));
        %figure;imagesc(f);axis image
        tmp = normxcorr2(f,I1);
        if blackAndWhite
            % Also detect black bumps
            f = imcomplement(f);
            tmp = max(tmp,normxcorr2(f,I1));
        end
        offset = floor((sz-1)/2);
        cdetect(:,:,i) = tmp(1+offset:end-offset,1+offset:end-offset);
    %     subaxis(2,3,i);
    %     imagesc(abs(cdetect(:,:,i))>0.7); axis image; colormap gray
    end
    [tmp_corrs,tmp_radius] = max(cdetect,[],3);
    cc = (tmp_corrs > thresh);
    %figure; imagesc(cc);axis image;colormap gray
    centroids = regionprops(bwlabel(cc,4),'Centroid');
    centroids = {centroids.Centroid};
    corrs = cellfun(@(x) tmp_corrs(uint16(x(2)),uint16(x(1))),centroids);
    radius = cellfun(@(x) sigma(tmp_radius(uint16(x(2)),uint16(x(1)))),centroids);
    centroids = reshape(cell2mat(centroids),2,size(centroids,2))';
end