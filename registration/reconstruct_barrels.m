function Irec = reconstruct_barrels1(I,M,debug)
    % This function aggregates the sections in I by mixing the gradients
    % of the sections and reconstructing the merged image.
    % We must mix images of the same nature. Thus the mixed gradient of
    % images i and i+1 cannot be mixed directly with image i+2.
    % 1) Mix pairs of images
    if debug
        figure
    end
    Irec1{1} = I{1};
    for i=1:length(I)-1
       Irec1{i} = mat2gray(mix_gradient(I{i},M{i},I{i+1},M{i+1},false));
       if debug
           subaxis(2,3,i);imagesc(Irec1{i});axis image;colormap gray
       end
    end
    % 2) Mix the reconstructed mixed gradient images obtained previously
    if debug
        figure
    end
    Irec2{1} = Irec1{1};
    for i=2:length(Irec1)
       Irec2{i} = mat2gray(mix_gradient(Irec2{i-1},M{i},Irec1{i},M{i+1},false));
       if debug
           subaxis(2,3,i-1);imagesc(Irec2{i});axis image;colormap gray
       end
    end
    Irec = Irec2{end};
end