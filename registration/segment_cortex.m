function [Iclean, BW] = segment_cortex(I,nClasses)
    [N,M] = size(I);
    centroid_idx = kmeans(double(I(:)),nClasses,'distance','sqEuclidean',...
                        'Replicates',3,'start','uniform');
    centroid_idx = reshape(centroid_idx,N,M);
    tmp = zeros(N,M);
    c = centroid_idx(floor(N/2),floor(M/2));
    tmp(centroid_idx == c) = 1;
    se = strel('disk', 5);
	tmp = imopen(tmp, se);
    % tmp = imclearborder(tmp);
    % tmp = imfill(tmp,'holes');
    Iseg = bwlabel(tmp);
    Ibin = zeros(N,M,'double');
    Ibin(Iseg == Iseg(floor(N/2),floor(M/2))) = 1;
    Ibin = imfill(Ibin,'holes');
    % Keep only the biggest connected component
    CC = bwconncomp(Ibin);
    S = regionprops(CC,'Centroid');
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    BW = zeros(N,M,'double');
    BW(CC.PixelIdxList{idx}) = 1;
    %% Mask input image
    Iclean = I.*BW;
end